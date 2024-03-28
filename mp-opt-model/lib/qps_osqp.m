function [x, f, eflag, output, lambda] = qps_osqp(H, c, A, l, u, xmin, xmax, x0, opt)
% qps_osqp - Quadratic Program Solver based on OSQP.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QPS_OSQP(H, C, A, L, U, XMIN, XMAX, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = QPS_OSQP(PROBLEM)
%   A wrapper function providing a standardized interface for using
%   OSQP to solve the following QP (quadratic programming) problem:
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
%               3 = even more verbose progress output
%           osqp_opt - options struct for OSQP (see
%               https://osqp.org/docs/interfaces/solver_settings.html),
%               value in verbose overrides these options
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : OSQP exit flag
%           1 = converged
%           0 or negative values OSQP status value
%           (see OSQP documentation for details)
%       OUTPUT : OSQP results struct
%           (see OSQP documentation for details)
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
%           qps_osqp(H, c, A, l, u, xmin, xmax, x0, opt)
%
%       x = qps_osqp(H, c, A, l, u)
%       x = qps_osqp(H, c, A, l, u, xmin, xmax)
%       x = qps_osqp(H, c, A, l, u, xmin, xmax, x0)
%       x = qps_osqp(H, c, A, l, u, xmin, xmax, x0, opt)
%       x = qps_osqp(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = qps_osqp(...)
%       [x, f] = qps_osqp(...)
%       [x, f, exitflag] = qps_osqp(...)
%       [x, f, exitflag, output] = qps_osqp(...)
%       [x, f, exitflag, output, lambda] = qps_osqp(...)
%
%
%   Example: (problem from from https://v8doc.sas.com/sashtml/iml/chap8/sect12.htm)
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
%       [x, f, s, out, lambda] = qps_osqp(H, c, A, l, u, xmin, [], x0, opt);
%
% See also qps_master, osqp.

%   MP-Opt-Model
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% check for OSQP
% if ~have_feature('osqp')
%     error('qps_osqp: requires OSQP (https://osqp.org)');
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
    lpqp = 'LP';
    if isempty(A) && isempty(xmin) && isempty(xmax)
        error('qps_osqp: LP problem must include constraints or variable bounds');
    else
        if ~isempty(A)
            nx = size(A, 2);
        elseif ~isempty(xmin)
            nx = length(xmin);
        else    % if ~isempty(xmax)
            nx = length(xmax);
        end
    end
else
    lpqp = 'QP';
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
if isempty(u)                   %% By default, linear inequalities are ...
    u = Inf(nA, 1);             %% ... unbounded above and ...
end
if isempty(l)
    l = -Inf(nA, 1);            %% ... unbounded below.
end

%% default options
if ~isempty(opt) && isfield(opt, 'verbose') && ~isempty(opt.verbose)
    verbose = opt.verbose;
else
    verbose = 0;
end

%% set up options struct for OSQP
osqp_opt = struct();
if ~isempty(opt) && isfield(opt, 'osqp_opt') && ~isempty(opt.osqp_opt)
    osqp_opt = nested_struct_copy(osqp_opt, opt.osqp_opt);
end
if verbose > 1
    osqp_opt.verbose = 1;
%     if verbose > 2
%         osqp_opt.
%     else
%         osqp_opt.
%     end
else
    osqp_opt.verbose = 0;
end

%% add variable bounds to linear constraints
if isempty(xmin)
    if isempty(xmax)
        k = [];
    else
        k = find(xmax < Inf);
    end
else
    if isempty(xmax)
        k = find(xmin > -Inf);
    else
        k = find(xmin > -Inf | xmax < Inf);
    end
end
nv = length(k);
if nv
    Av = sparse(1:nv, k, 1, nv, nx);
    if isempty(xmin)
        lv = -Inf(nv,1);
    else
        lv = xmin(k);
    end
    if isempty(xmax)
        uv = Inf(nv,1);
    else
        uv = xmax(k);
    end
    AA = [A; Av];
    ll = [l; lv];
    uu = [u; uv];
else
    AA = A;
    ll = l;
    uu = u;
end

%% set up model
o = osqp();
o.setup(H, c, AA, ll, uu, osqp_opt);
if ~isempty(x0)
    o.warm_start('x', x0);
end

%% solve model
if verbose
    vn = osqpver;
    fprintf('OSQP Version %s -- %s solver\n', vn, lpqp);
end
res = o.solve();
if verbose
    fprintf('OSQP solution status: %s\n', res.info.status);
end

%% extract results
x = res.x;
f = res.info.obj_val;
eflag = res.info.status_val;
if eflag > 1
    eflag = -10 * eflag;
end
output = res;

%% separate duals into binding lower & upper bounds
mu_l = -res.y;
mu_u =  res.y;
mu_l(mu_l < 0) = 0;
mu_u(mu_u < 0) = 0;

%% variable bounds
lam_lower = zeros(nx, 1);
lam_lower(k) = mu_l(nA+1:end);
lam_upper = zeros(nx, 1);
lam_upper(k) = mu_u(nA+1:end);

%% package lambdas
lambda = struct( ...
    'mu_l', mu_l(1:nA), ...
    'mu_u', mu_u(1:nA), ...
    'lower', lam_lower, ...
    'upper', lam_upper ...
);
