function [x, f, eflag, output, lambda] = miqps_gurobi(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
%MIQPS_GUROBI  Mixed Integer Quadratic Program Solver based on GUROBI.
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       MIQPS_GUROBI(H, C, A, L, U, XMIN, XMAX, X0, OPT)
%   A wrapper function providing a MATPOWER standardized interface for using
%   GUROBI to solve the following QP (quadratic programming)
%   problem:
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
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%               3 = even more verbose progress output
%           skip_prices (0) - flag that specifies whether or not to
%               skip the price computation stage, in which the problem
%               is re-solved for only the continuous variables, with all
%               others being constrained to their solved values
%           price_stage_warn_tol (1e-7) - tolerance on the objective fcn
%               value and primal variable relative match required to avoid
%               mis-match warning message
%           grb_opt - options struct for GUROBI, value in
%               verbose overrides these options
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, vtype, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : GUROBI exit flag
%           1 = converged
%           0 or negative values = negative of GUROBI exit flag
%           (see GUROBI documentation for details)
%       OUTPUT : GUROBI output struct
%           (see GUROBI documentation for details)
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
%           miqps_gurobi(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
%
%       x = miqps_gurobi(H, c, A, l, u)
%       x = miqps_gurobi(H, c, A, l, u, xmin, xmax)
%       x = miqps_gurobi(H, c, A, l, u, xmin, xmax, x0)
%       x = miqps_gurobi(H, c, A, l, u, xmin, xmax, x0, vtype)
%       x = miqps_gurobi(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
%       x = miqps_gurobi(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, vtype, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = miqps_gurobi(...)
%       [x, f] = miqps_gurobi(...)
%       [x, f, exitflag] = miqps_gurobi(...)
%       [x, f, exitflag, output] = miqps_gurobi(...)
%       [x, f, exitflag, output, lambda] = miqps_gurobi(...)
%
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
%       [x, f, s, out, lambda] = miqps_gurobi(H, c, A, l, u, xmin, [], x0, vtype, opt);
%
%   See also GUROBI.

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

%% define nx, set default values for missing optional inputs
if isempty(H) || ~any(any(H))
    if isempty(A) && isempty(xmin) && isempty(xmax)
        error('miqps_gurobi: LP problem must include constraints or variable bounds');
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
if isempty(xmin)                %% By default, optimization variables are ...
    xmin = -Inf(nx, 1);         %% ... unbounded below and ...
end
if isempty(xmax)
    xmax = Inf(nx, 1);          %% ... unbounded above.
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

%% set up options struct for Gurobi
if ~isempty(opt) && isfield(opt, 'grb_opt') && ~isempty(opt.grb_opt)
    g_opt = gurobi_options(opt.grb_opt);
else
    g_opt = gurobi_options;
end
if verbose > 1
    g_opt.LogToConsole = 1;
    g_opt.OutputFlag = 1;
    if verbose > 2
        g_opt.DisplayInterval = 1;
    else
        g_opt.DisplayInterval = 100;
    end
else
    g_opt.LogToConsole = 0;
    g_opt.OutputFlag = 0;
end

if ~issparse(A)
    A = sparse(A);
end
if issparse(c);
    c = full(c);
end

%% split up linear constraints
ieq = find( abs(u-l) <= eps );          %% equality
igt = find( u >=  1e10 & l > -1e10 );   %% greater than, unbounded above
ilt = find( l <= -1e10 & u <  1e10 );   %% less than, unbounded below
ibx = find( (abs(u-l) > eps) & (u < 1e10) & (l > -1e10) );

%% grab some dimensions
nlt = length(ilt);      %% number of upper bounded linear inequalities
ngt = length(igt);      %% number of lower bounded linear inequalities
nbx = length(ibx);      %% number of doubly bounded linear inequalities
neq = length(ieq);      %% number of equalities
niq = nlt+ngt+2*nbx;    %% number of inequalities

%% set up model
m.A     = [ A(ieq, :); A(ilt, :); -A(igt, :); A(ibx, :); -A(ibx, :) ];
m.rhs   = [ u(ieq);    u(ilt);    -l(igt);    u(ibx);    -l(ibx)    ];
m.sense = char([ double('=')*ones(1,neq) double('<')*ones(1,niq) ]);
m.lb = xmin;
m.ub = xmax;
m.obj = c';
if ~isempty(vtype)
    m.vtype = vtype;
end
if isempty(vtype) || isempty(find(vtype == 'B' | vtype == 'I' | ...
        vtype == 'S' | vtype == 'N'))
    mi = 0;
else
    mi = 1;
end

%% call the solver
if isempty(H) || ~any(any(H))
    lpqp = 'LP';
else
    lpqp = 'QP';
    if ~issparse(H)
        H = sparse(H);
    end
    m.Q = 0.5 * H;
end
if mi
    lpqp = ['MI' lpqp];
end
if verbose
    alg_names = {
        'automatic',
        'primal simplex',
        'dual simplex',
        'interior point',
        'concurrent',
        'deterministic concurrent'
    };
    vn = gurobiver;
    fprintf('Gurobi Version %s -- %s %s solver\n', ...
        vn, alg_names{g_opt.Method+2}, lpqp);
end
results = gurobi(m, g_opt);
switch results.status
    case 'LOADED',          %% 1
        eflag = -1;
    case 'OPTIMAL',         %% 2, optimal solution found
        eflag = 1;
    case 'INFEASIBLE',      %% 3
        eflag = -3;
    case 'INF_OR_UNBD',     %% 4
        eflag = -4;
    case 'UNBOUNDED',       %% 5
        eflag = -5;
    case 'CUTOFF',          %% 6
        eflag = -6;
    case 'ITERATION_LIMIT', %% 7
        eflag = -7;
    case 'NODE_LIMIT',      %% 8
        eflag = -8;
    case 'TIME_LIMIT',      %% 9
        eflag = -9;
    case 'SOLUTION_LIMIT',  %% 10
        eflag = -10;
    case 'INTERRUPTED',     %% 11
        eflag = -11;
    case 'NUMERIC',         %% 12
        eflag = -12;
    case 'SUBOPTIMAL',      %% 13
        eflag = -13;
    case 'INPROGRESS',      %% 14
        eflag = -14;
    otherwise,
        eflag = 0;
end
output = results;

%% check for empty results (in case optimization failed)
if ~isfield(results, 'x') || isempty(results.x)
    x = NaN(nx, 1);
    lam.lower   = NaN(nx, 1);
    lam.upper   = NaN(nx, 1);
else
    x = results.x;
    lam.lower   = zeros(nx, 1);
    lam.upper   = zeros(nx, 1);
end
if ~isfield(results, 'objval') || isempty(results.objval)
    f = NaN;
else
    f = results.objval;
end
if ~isfield(results, 'pi') || isempty(results.pi)
    pi  = NaN(length(m.rhs), 1);
else
    pi  = results.pi;
end
if ~isfield(results, 'rc') || isempty(results.rc)
    rc  = NaN(nx, 1);
else
    rc  = results.rc;
end

kl = find(rc > 0);   %% lower bound binding
ku = find(rc < 0);   %% upper bound binding
lam.lower(kl)   =  rc(kl);
lam.upper(ku)   = -rc(ku);
lam.eqlin   = pi(1:neq);
lam.ineqlin = pi(neq+(1:niq));
mu_l        = zeros(nA, 1);
mu_u        = zeros(nA, 1);

%% repackage lambdas
kl = find(lam.eqlin > 0);   %% lower bound binding
ku = find(lam.eqlin < 0);   %% upper bound binding

mu_l(ieq(kl)) = lam.eqlin(kl);
mu_l(igt) = -lam.ineqlin(nlt+(1:ngt));
mu_l(ibx) = -lam.ineqlin(nlt+ngt+nbx+(1:nbx));

mu_u(ieq(ku)) = -lam.eqlin(ku);
mu_u(ilt) = -lam.ineqlin(1:nlt);
mu_u(ibx) = -lam.ineqlin(nlt+ngt+(1:nbx));

lambda = struct( ...
    'mu_l', mu_l, ...
    'mu_u', mu_u, ...
    'lower', lam.lower, ...
    'upper', lam.upper ...
);

if mi && eflag == 1 && (~isfield(opt, 'skip_prices') || ~opt.skip_prices)
    if verbose
        fprintf('--- Integer stage complete, starting price computation stage ---\n');
    end
    if isfield(opt, 'price_stage_warn_tol') && ~isempty(opt.price_stage_warn_tol)
        tol = opt.price_stage_warn_tol;
    else
        tol = 1e-7;
    end
    if length(vtype) == 1
        if vtype == 'I' || vtype == 'B' || vtype == 'N'
            k = (1:nx);
        elseif vtype == 'S'
            k = find(x == 0);
        end
    else    %% length(vtype) == nx
        k = find(vtype == 'I' | vtype == 'B' | vtype == 'N' | ...
                (vtype == 'S' & x' == 0));
    end
    
    x(k) = round(x(k));
    xmin(k) = x(k);
    xmax(k) = x(k);
    x0 = x;
%     opt.grb_opt.Method = 0;     %% primal simplex
    
    [x_, f_, eflag_, output_, lambda] = qps_gurobi(H, c, A, l, u, xmin, xmax, x0, opt);
    if eflag ~= eflag_
        error('miqps_gurobi: EXITFLAG from price computation stage = %d', eflag_);
    end
    if abs(f - f_)/max(abs(f), 1) > tol
        warning('miqps_gurobi: relative mismatch in objective function value from price computation stage = %g', abs(f - f_)/max(abs(f), 1));
    end
    xn = x;
    xn(abs(xn)<1) = 1;
    [mx, k] = max(abs(x - x_) ./ xn);
    if mx > tol
        warning('miqps_gurobi: max relative mismatch in x from price computation stage = %g (%g)', mx, x(k));
    end
    output.price_stage = output_;
end
