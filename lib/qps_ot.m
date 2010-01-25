function [x, f, eflag, output, lambda] = qps_ot(H, c, A, l, u, xmin, xmax, x0, opt)
%QPS_OT  Quadratic Program Solver based on quadprog()/linprog()
%   A wrapper function providing a MATPOWER standardized interface for using
%   quadprog() or linprog() from the Optimization Toolbox to solve the
%   following QP (quadratic programming) problem:
%
%       min 1/2 x'*H*x + c'*x
%        x
%
%   subject to
%
%       l <= A*x <= u       (linear constraints)
%       xmin <= x <= xmax   (variable bounds)
%
%   [x, fval, exitflag, output, lambda] = ...
%       qps_ot(H, c, A, l, u, xmin, xmax, x0, opt)
%
%   x = qps_ot(H, c, A, l, u)
%   x = qps_ot(H, c, A, l, u, xmin, xmax)
%   x = qps_ot(H, c, A, l, u, xmin, xmax, x0)
%   x = qps_ot(H, c, A, l, u, xmin, xmax, x0, opt)
%   x = qps_ot(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'f' and 'x0' are optional
%   x = qps_ot(...)
%   [x, fval] = qps_ot(...)
%   [x, fval, exitflag] = qps_ot(...)
%   [x, fval, exitflag, output] = qps_ot(...)
%   [x, fval, exitflag, output, lambda] = qps_ot(...)
%
%   Inputs:
%       H : matrix (possibly sparse) of quadratic cost coefficients
%       c : vector of linear cost coefficients
%       A, l, u : define the optional linear constraints. Default
%           values for the elements of l and u are -Inf and Inf,
%           respectively.
%       xmin, xmax : optional lower and upper bounds on the
%           x variables, defaults are -Inf and Inf, respectively.
%       x0 : optional starting value of optimization vector x
%       opt : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           max_it (0) - maximum number of iterations allowed
%               0 = use algorithm default
%           ot_opt - options struct for quadprog()/linprog(), values in
%               verbose and max_it override these options
%       problem : The inputs can alternatively be supplied in a single
%           struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, opt
%
%   Outputs:
%       x : solution vector
%       fval : final objective function value
%       exitflag : quadprog()/linprog() exit flag
%           (see quadprog and linprog documentation for details)
%       output : quadprog()/linprog() output structure
%           (see quadprog and linprog documentation for details)
%       lambda : struct containing the Langrange and Kuhn-Tucker
%           multipliers on the constraints, with fields:
%           mu_l - lower bound on linear constraints
%           mu_u - upper bound on linear constraints
%           lower - lower bound on optimization variables
%           upper - upper bound on optimization variables
%
%   Note the calling syntax is almost identical to that of 'quadprog'
%   from MathWorks' Optimization Toolbox. The main difference is that
%   the linear constraints are specified with A, l, u instead of
%   A, b, Aeq, beq.
%
%   See also quadprog and linprog

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2010 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% check for Optimization Toolbox
% if ~have_fcn('quadprog')
%     error('qps_ot: requires the Optimization Toolbox');
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
        error('qps_ot: LP problem must include constraints or variable bounds');
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
if  ~isempty(A) && (isempty(l) || all(l == -Inf)) && ...
                   (isempty(u) || all(u == Inf))
    A = sparse(0,nx);           %% no limits => no linear constraints
end
nA = size(A, 1);                %% number of original linear constraints
if isempty(u)                   %% By default, linear inequalities are ...
    u = Inf * ones(nA, 1);      %% ... unbounded above and ...
end
if isempty(l)
    l = -Inf * ones(nA, 1);     %% ... unbounded below.
end
if isempty(xmin)                %% By default, optimization variables are ...
    xmin = -Inf * ones(nx, 1);  %% ... unbounded below and ...
end
if isempty(xmax)
    xmax = Inf * ones(nx, 1);   %% ... unbounded above.
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
if ~isempty(opt) && isfield(opt, 'max_it') && ~isempty(opt.max_it)
    max_it = opt.max_it;
else
    max_it = 0;
end

%% split up linear constraints
ieq = find( abs(u-l) <= eps );          %% equality
igt = find( u >=  1e10 & l > -1e10 );   %% greater than, unbounded above
ilt = find( l <= -1e10 & u <  1e10 );   %% less than, unbounded below
ibx = find( (abs(u-l) > eps) & (u < 1e10) & (l > -1e10) );
Ae = A(ieq, :);
be = u(ieq);
Ai  = [ A(ilt, :); -A(igt, :); A(ibx, :); -A(ibx, :) ];
bi  = [ u(ilt);    -l(igt);    u(ibx);    -l(ibx)];

%% grab some dimensions
nlt = length(ilt);      %% number of upper bounded linear inequalities
ngt = length(igt);      %% number of lower bounded linear inequalities
nbx = length(ibx);      %% number of doubly bounded linear inequalities

%% set up options
if ~isempty(opt) && isfield(opt, 'ot_opt') && ~isempty(opt.ot_opt)
    ot_opt = opt.ot_opt;
else
    if isempty(H) || ~any(any(H))
        ot_opt = optimset('linprog');
    else
        ot_opt = optimset('quadprog');
    end
end
if max_it
    ot_opt = optimset(ot_opt, 'MaxIter', max_it);
end
if verbose > 1
    ot_opt = optimset(ot_opt, 'Display', 'iter');   %% seems to be same as 'final'
elseif verbose == 1
    ot_opt = optimset(ot_opt, 'Display', 'final');
else
    ot_opt = optimset(ot_opt, 'Display', 'off');
end

%% call the solver
if isempty(H) || ~any(any(H))
    [x, f, eflag, output, lam] = ...
        linprog(c, Ai, bi, Ae, be, xmin, xmax, x0, ot_opt);
else
    [x, f, eflag, output, lam] = ...
        quadprog(H, c, Ai, bi, Ae, be, xmin, xmax, x0, ot_opt);
end

%% repackage lambdas
kl = find(lam.eqlin < 0);   %% lower bound binding
ku = find(lam.eqlin > 0);   %% upper bound binding

mu_l = zeros(nA, 1);
mu_l(ieq(kl)) = -lam.eqlin(kl);
mu_l(igt) = lam.ineqlin(nlt+(1:ngt));
mu_l(ibx) = lam.ineqlin(nlt+ngt+nbx+(1:nbx));

mu_u = zeros(nA, 1);
mu_u(ieq(ku)) = lam.eqlin(ku);
mu_u(ilt) = lam.ineqlin(1:nlt);
mu_u(ibx) = lam.ineqlin(nlt+ngt+(1:nbx));

lambda = struct( ...
    'mu_l', mu_l, ...
    'mu_u', mu_u, ...
    'lower', lam.lower, ...
    'upper', lam.upper ...
);
