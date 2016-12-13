function [x, f, eflag, output, lambda] = qps_ipopt(H, c, A, l, u, xmin, xmax, x0, opt)
%QPS_IPOPT  Quadratic Program Solver based on IPOPT.
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QPS_IPOPT(H, C, A, L, U, XMIN, XMAX, X0, OPT)
%   Uses IPOPT to solve the following QP (quadratic programming) problem:
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
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           ipopt_opt - options struct for IPOPT, value in verbose
%               overrides these options
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : exit flag
%           1 = first order optimality conditions satisfied
%           0 = maximum number of iterations reached
%           -1 = numerically failed
%       OUTPUT : output struct with the following fields:
%           iterations - number of iterations performed
%           hist - struct array with trajectories of the following:
%                   feascond, gradcond, compcond, costcond, gamma,
%                   stepsize, obj, alphap, alphad
%           message - exit message
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
%           qps_ipopt(H, c, A, l, u, xmin, xmax, x0, opt)
%
%       x = qps_ipopt(H, c, A, l, u)
%       x = qps_ipopt(H, c, A, l, u, xmin, xmax)
%       x = qps_ipopt(H, c, A, l, u, xmin, xmax, x0)
%       x = qps_ipopt(H, c, A, l, u, xmin, xmax, x0, opt)
%       x = qps_ipopt(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = qps_ipopt(...)
%       [x, f] = qps_ipopt(...)
%       [x, f, exitflag] = qps_ipopt(...)
%       [x, f, exitflag, output] = qps_ipopt(...)
%       [x, f, exitflag, output, lambda] = qps_ipopt(...)
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
%       [x, f, s, out, lambda] = qps_ipopt(H, c, A, l, u, xmin, [], x0, opt);
%
%   See also IPOPT, IPOPT_OPTIONS.
%   http://www.coin-or.org/projects/Ipopt.xml.

%   MATPOWER
%   Copyright (c) 2010-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% check for IPOPT
% if ~have_fcn('ipopt')
%     error('qps_ipopt: requires IPOPT (http://www.coin-or.org/projects/Ipopt.xml)');
% end

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

%% define nx, set default values for missing optional inputs
if isempty(H) || ~any(any(H))
    if isempty(A) && isempty(xmin) && isempty(xmax)
        error('qps_ipopt: LP problem must include constraints or variable bounds');
    else
        if ~isempty(A)
            nx = size(A, 2);
        elseif ~isempty(xmin)
            nx = length(xmin);
        else    % if ~isempty(xmax)
            nx = length(xmax);
        end
    end
    H = sparse(nx,nx);
else
    nx = size(H, 1);
end
if isempty(c)
    c = zeros(nx, 1);
end
if isempty(A) || (~isempty(A) && (isempty(l) || all(l == -Inf)) && ...
                                 (isempty(u) || all(u == Inf)))
    A = sparse(0,nx);           %% no limits => no linear constraints
end
nA = size(A, 1);                %% number of original linear constraints
if nA
    if isempty(u)               %% By default, linear inequalities are ...
        u = Inf(nA, 1);             %% ... unbounded above and ...
    end
    if isempty(l)
        l = -Inf(nA, 1);            %% ... unbounded below.
    end
end
if isempty(x0)
    x0 = zeros(nx, 1);
end

%% default options
if ~isempty(opt) && isfield(opt, 'verbose') && ~isempty(opt.verbose)
    verbose = opt.verbose;
else
    verbose = 0;
end

%% make sure args are sparse/full as expected by IPOPT
if ~isempty(H)
    if ~issparse(H)
        H = sparse(H);
    end
end
if ~issparse(A)
    A = sparse(A);
end


%%-----  run optimization  -----
%% set options struct for IPOPT
if ~isempty(opt) && isfield(opt, 'ipopt_opt') && ~isempty(opt.ipopt_opt)
    options.ipopt = ipopt_options(opt.ipopt_opt);
else
    options.ipopt = ipopt_options;
end
options.ipopt.jac_c_constant    = 'yes';
options.ipopt.jac_d_constant    = 'yes';
options.ipopt.hessian_constant  = 'yes';
options.ipopt.least_square_init_primal  = 'yes';
options.ipopt.least_square_init_duals   = 'yes';
% options.ipopt.mehrotra_algorithm        = 'yes';        %% default 'no'
if verbose
    options.ipopt.print_level = min(12, verbose*2+1);
else
    options.ipopt.print_level = 0;
end

%% define variable and constraint bounds, if given
if nA
    options.cu = u;
    options.cl = l;
end
if ~isempty(xmin)
    options.lb = xmin;
end
if ~isempty(xmax)
    options.ub = xmax;
end

%% assign function handles
funcs.objective         = @(x) 0.5 * x' * H * x + c' * x;
funcs.gradient          = @(x) H * x + c;
funcs.constraints       = @(x) A * x;
funcs.jacobian          = @(x) A;
funcs.jacobianstructure = @() A;
funcs.hessian           = @(x, sigma, lambda) tril(H);
funcs.hessianstructure  = @() tril(H);

%% run the optimization
[x, info] = ipopt(x0,funcs,options);

if info.status == 0 || info.status == 1
    eflag = 1;
else
    eflag = 0;
end
if isfield(info, 'iter')
    output.iterations = info.iter;
end
output.info       = info.status;
f = funcs.objective(x);

%% check for empty results (in case optimization failed)
if isempty(info.lambda)
	lam = NaN(nA, 1);
else
	lam = info.lambda;
end
if isempty(info.zl) && ~isempty(xmin)
	zl = NaN(nx, 1);
else
	zl = info.zl;
end
if isempty(info.zu) && ~isempty(xmax)
	zu = NaN(nx, 1);
else
	zu = info.zu;
end

%% repackage lambdas
kl = find(lam < 0);                     %% lower bound binding
ku = find(lam > 0);                     %% upper bound binding
mu_l = zeros(nA, 1);
mu_l(kl) = -lam(kl);
mu_u = zeros(nA, 1);
mu_u(ku) = lam(ku);

lambda = struct( ...
    'mu_l', mu_l, ...
    'mu_u', mu_u, ...
    'lower', zl, ...
    'upper', zu    );
