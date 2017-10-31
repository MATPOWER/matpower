function [x, f, eflag, output, lambda] = miqps_matpower(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
%MIQPS_MATPOWER  Mixed Integer Quadratic Program Solver for MATPOWER.
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       MIQPS_MATPOWER(H, C, A, L, U, XMIN, XMAX, X0, VTYPE, OPT)
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
%       VTYPE : character string of length NX (number of elements in X),
%               or 1 (value applies to all variables in x),
%               allowed values are 'C' (continuous), 'B' (binary),
%               'I' (integer), 'S' (semi-continuous), or 'N' (semi-integer).
%               (MOSEK, GLPK, OT allow only 'C', 'B', or 'I')
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           alg ('DEFAULT') : determines which solver to use, can be either
%                   a string (new-style) or a numerical alg code (old-style)
%               'DEFAULT' : (or 0) automatic, first available of CPLEX,
%                       Gurobi, MOSEK, Opt Tbx (MILPs only), GLPK (MILPs only)
%               'CPLEX'   : (or 500) CPLEX
%               'GLPK'    : GLPK, (MILP problems only, i.e. empty H matrix)
%               'GUROBI'  : (or 700) Gurobi
%               'MOSEK'   : (or 600) MOSEK
%               'OT'      : (or 300) Optimization Toolbox, INTLINPROG
%                           (MILP problems only, i.e. empty H matrix)
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           skip_prices (0) - flag that specifies whether or not to
%               skip the price computation stage, in which the problem
%               is re-solved for only the continuous variables, with all
%               others being constrained to their solved values
%           price_stage_warn_tol (1e-7) - tolerance on the objective fcn
%               value and primal variable relative match required to avoid
%               mis-match warning message
%           cplex_opt - options struct for CPLEX
%           glpk_opt    - options struct for GLPK
%           grb_opt   - options struct for GBUROBI_MEX
%           intlinprog_opt - options struct for INTLINPROG
%           linprog_opt - options struct for LINPROG
%           mosek_opt - options struct for MOSEK
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, vtype, opt
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
%           miqps_matpower(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
%
%       x = miqps_matpower(H, c, A, l, u)
%       x = miqps_matpower(H, c, A, l, u, xmin, xmax)
%       x = miqps_matpower(H, c, A, l, u, xmin, xmax, x0)
%       x = miqps_matpower(H, c, A, l, u, xmin, xmax, x0, vtype)
%       x = miqps_matpower(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
%       x = miqps_matpower(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, vtype, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = miqps_matpower(...)
%       [x, f] = miqps_matpower(...)
%       [x, f, exitflag] = miqps_matpower(...)
%       [x, f, exitflag, output] = miqps_matpower(...)
%       [x, f, exitflag, output, lambda] = miqps_matpower(...)
%
%   Example: (problem from from http://www.jmu.edu/docs/sasdoc/sashtml/iml/chap8/sect12.htm)
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
%       [x, f, s, out, lambda] = miqps_matpower(H, c, A, l, u, xmin, [], x0, vtype, opt);

%   MATPOWER
%   Copyright (c) 2010-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- input argument handling  -----
%% gather inputs
if nargin == 1 && isstruct(H)       %% problem struct
    p = H;
    if isfield(p, 'opt'),   opt = p.opt;    else,   opt = [];   end
    if isfield(p, 'vtype'), vtype = p.vtype;else,   vtype = []; end
    if isfield(p, 'x0'),    x0 = p.x0;      else,   x0 = [];    end
    if isfield(p, 'xmax'),  xmax = p.xmax;  else,   xmax = [];  end
    if isfield(p, 'xmin'),  xmin = p.xmin;  else,   xmin = [];  end
    if isfield(p, 'u'),     u = p.u;        else,   u = [];     end
    if isfield(p, 'l'),     l = p.l;        else,   l = [];     end
    if isfield(p, 'A'),     A = p.A;        else,   A = [];     end
    if isfield(p, 'c'),     c = p.c;        else,   c = [];     end
    if isfield(p, 'H'),     H = p.H;        else,   H = [];     end
else                                %% individual args
    if nargin < 10
        opt = [];
        if nargin < 9
            vtype = [];
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
end

%% default options
if ~isempty(opt) && isfield(opt, 'alg') && ~isempty(opt.alg)
    alg = opt.alg;
    %% convert integer codes to string values
    if ~ischar(alg)
        switch alg
            case 0
                alg = 'DEFAULT';
            case 300
                alg = 'OT';
            case 500
                alg = 'CPLEX';
            case 600
                alg = 'MOSEK';
            case 700
                alg = 'GUROBI';
            otherwise
                error('miqps_matpower: %d is not a valid algorithm code', alg);
        end
    end
else
    alg = 'DEFAULT';
end
if ~isempty(opt) && isfield(opt, 'verbose') && ~isempty(opt.verbose)
    verbose = opt.verbose;
else
    verbose = 0;
end
if strcmp(alg, 'DEFAULT')
    if have_fcn('gurobi')       %% use Gurobi by default, if available
        alg = 'GUROBI';
    elseif have_fcn('cplex')    %% if not, then CPLEX, if available
        alg = 'CPLEX';
    elseif have_fcn('mosek')    %% if not, then MOSEK, if available
        alg = 'MOSEK';
    elseif isempty(H) || ~any(any(H))   %% if not, and linear objective
        if have_fcn('intlinprog')       %% then Optimization Tbx, if available
            alg = 'OT';
        elseif have_fcn('glpk')         %% if not, and then GLPK, if available
            alg = 'GLPK';
        end
    else
        error('miqps_matpower: no solvers available - requires CPLEX, Gurobi, MOSEK, INTLINPROG or GLPK');
    end
end

%%----- call the appropriate solver  -----
switch alg
    case 'CPLEX'
        [x, f, eflag, output, lambda] = ...
            miqps_cplex(H, c, A, l, u, xmin, xmax, x0, vtype, opt);
    case 'GLPK'
        [x, f, eflag, output, lambda] = ...
            miqps_glpk(H, c, A, l, u, xmin, xmax, x0, vtype, opt);
    case 'GUROBI'
        [x, f, eflag, output, lambda] = ...
            miqps_gurobi(H, c, A, l, u, xmin, xmax, x0, vtype, opt);
    case 'MOSEK'
        [x, f, eflag, output, lambda] = ...
            miqps_mosek(H, c, A, l, u, xmin, xmax, x0, vtype, opt);
    case 'OT'
        [x, f, eflag, output, lambda] = ...
            miqps_ot(H, c, A, l, u, xmin, xmax, x0, vtype, opt);
    otherwise
        error('miqps_matpower: ''%s'' is not a valid algorithm code', alg);
end
if ~isfield(output, 'alg') || isempty(output.alg)
    output.alg = alg;
end
