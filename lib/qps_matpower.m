function [x, f, eflag, output, lambda] = qps_matpower(H, c, A, l, u, xmin, xmax, x0, opt)
%QPS_MATPOWER  Quadratic Program Solver for MATPOWER.
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QPS_MATPOWER(H, C, A, L, U, XMIN, XMAX, X0, OPT)
%   A common wrapper function for various QP solvers. 
%   Solves the following QP (quadratic programming) problem:
%
%       min 1/2 X'*H*X + C'*X
%        X
%
%   subject to
%
%       L <= A*X <= U       (linear constraints)
%       XMIN <= X <= XMAX   (variable bounds)
%
%   Inputs (all optional except H, C, A and L):
%       H : matrix (possibly sparse) of quadratic cost coefficients
%       C : vector of linear cost coefficients
%       A, L, U : define the optional linear constraints. Default
%           values for the elements of L and U are -Inf and Inf,
%           respectively.
%       XMIN, XMAX : optional lower and upper bounds on the
%           X variables, defaults are -Inf and Inf, respectively.
%       X0 : optional starting value of optimization vector X
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           alg (0) - determines which solver to use
%                 0 = use first available solver
%               100 = BPMPD_MEX
%               200 = MIPS, Matlab Interior Point Solver
%                     pure Matlab implementation of a primal-dual
%                     interior point method
%               250 = MIPS-sc, a step controlled variant of MIPS
%               300 = Optimization Toolbox, QUADPROG or LINPROG
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           max_it (0) - maximum number of iterations allowed
%               0 = use algorithm default
%           bp_opt - options vector for BP
%           mips_opt - options struct for QPS_MIPS
%           ot_opt - options struct for QUADPROG/LINPROG
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : exit flag
%           1 = converged
%           0 or negative values = algorithm specific failure codes
%       OUTPUT : output struct with the following fields:
%           alg - algorithm code of solver used
%           (others) - algorithm specific fields
%       LAMBDA : struct containing the Langrange and Kuhn-Tucker
%           multipliers on the constraints, with fields:
%           mu_l - lower (left-hand) limit on linear constraints
%           mu_u - upper (right-hand) limit on linear constraints
%           lower - lower bound on optimization variables
%           upper - upper bound on optimization variables
%
%   Note the calling syntax is almost identical to that of QUADPROG
%   from MathWorks' Optimization Toolbox. The main difference is that
%   the linear constraints are specified with A, L, U instead of
%   A, B, Aeq, Beq.
%
%   Calling syntax options:
%       [x, f, exitflag, output, lambda] = ...
%           qps_matpower(H, c, A, l, u, xmin, xmax, x0, opt)
%
%       x = qps_matpower(H, c, A, l, u)
%       x = qps_matpower(H, c, A, l, u, xmin, xmax)
%       x = qps_matpower(H, c, A, l, u, xmin, xmax, x0)
%       x = qps_matpower(H, c, A, l, u, xmin, xmax, x0, opt)
%       x = qps_matpower(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'H', 'c', 'A' and 'l' are optional
%       x = qps_matpower(...)
%       [x, f] = qps_matpower(...)
%       [x, f, exitflag] = qps_matpower(...)
%       [x, f, exitflag, output] = qps_matpower(...)
%       [x, f, exitflag, output, lambda] = qps_matpower(...)
%
%   Example: (problem from from http://www.uc.edu/sashtml/iml/chap8/sect12.htm)
%       H = [   1003.1  4.3     6.3     5.9;
%               4.3     2.2     2.1     3.9;
%               6.3     2.1     3.5     4.8;
%               5.9     3.9     4.8     10  ];
%       c = zeros(4,1);
%       A = [   1       1       1       1;
%               0.17    0.11    0.10    0.18    ];
%       l = [1; 0.10];
%       u = [1; Inf];
%       xmin = zeros(4,1);
%       x0 = [1; 0; 0; 1];
%       opt = struct('verbose', 2);
%       [x, f, s, out, lam] = qps_matpower(H, c, A, l, u, xmin, [], x0, opt);

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2010 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- input argument handling  -----
%% gather inputs
if nargin == 1 && isstruct(H)       %% problem struct
    p = H;
    if isfield(p, 'opt'),   opt = p.opt;    else,   opt = [];   end
    if isfield(p, 'x0'),    x0 = p.x0;      else,   x0 = [];    end
    if isfield(p, 'xmax'),  xmax = p.xmax;  else,   xmax = [];  end
    if isfield(p, 'xmin'),  xmin = p.xmin;  else,   xmin = [];  end
    if isfield(p, 'u'),     u = p.u;        else,   u = [];     end
    if isfield(p, 'l'),     l = p.l;        else,   l = [];     end
    if isfield(p, 'A'),     A = p.A;        else,   A = [];     end
    if isfield(p, 'c'),     c = p.c;        else,   c = [];     end
    if isfield(p, 'H'),     H = p.H;        else,   H = [];     end
else                                %% individual args
    if nargin < 9
        opt = [];
        if nargin < 8
            x0 = [];
            if nargin < 7
                xmax = [];
                if nargin < 6
                    xmin = [];
                end
            end
        end
    end
end

%% default options
if ~isempty(opt) && isfield(opt, 'alg') && ~isempty(opt.alg)
    alg = opt.alg;
else
    alg = 0;
end
if ~isempty(opt) && isfield(opt, 'verbose') && ~isempty(opt.verbose)
    verbose = opt.verbose;
else
    verbose = 0;
end
if alg == 0
    if have_fcn('bpmpd')    %% use BPMPD_MEX by default if available
        alg = 100;          %% BPMPD_MEX
    else
        alg = 200;          %% MIPS
    end
end

%%----- call the appropriate solver  -----
if alg == 100                       %% use BPMPD_MEX
    [x, f, eflag, output, lambda] = ...
        qps_bpmpd(H, c, A, l, u, xmin, xmax, x0, opt);

    if eflag == -99
        if verbose
            fprintf('         Retrying with QPS_MIPS solver ...\n\n');
        end
        %% save (incorrect) solution from BPMPD
        bpmpd = struct('x', x, 'f', f, 'eflag', eflag, ...
                        'output', output, 'lambda', lambda);
        opt.alg = 200;
        [x, f, eflag, output, lambda] = ...
            qps_matpower(H, c, A, l, u, xmin, xmax, x0, opt);
        output.bpmpd = bpmpd;
    end
elseif alg == 200 || alg == 250     %% use MIPS or sc-MIPS
    %% set up options
    if ~isempty(opt) && isfield(opt, 'mips_opt') && ~isempty(opt.mips_opt)
        mips_opt = opt.mips_opt;
    else
        mips_opt = [];
    end
    if ~isempty(opt) && isfield(opt, 'max_it') && ~isempty(opt.max_it)
        mips_opt.max_it = opt.max_it;
    end
    if alg == 200
        mips_opt.step_control = 0;
    else
        mips_opt.step_control = 1;
    end
    mips_opt.verbose = verbose;
    
    mlver = ver('matlab');
    if str2double(mlver.Version(1)) < 7    %% anonymous functions not available
        solver = 'qps_mips6';
    else
        solver = 'qps_mips';
    end

    %% call solver
    [x, f, eflag, output, lambda] = ...
        feval(solver, H, c, A, l, u, xmin, xmax, x0, mips_opt);
elseif alg == 300                   %% use QUADPROG or LINPROG from Opt Tbx ver 2.x+
    [x, f, eflag, output, lambda] = ...
        qps_ot(H, c, A, l, u, xmin, xmax, x0, opt);
else
    error('qps_matpower: %d is not a valid algorithm code', alg);
end
if ~isfield(output, 'alg') || isempty(output.alg)
    output.alg = alg;
end
