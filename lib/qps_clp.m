function [x, f, eflag, output, lambda] = qps_clp(H, c, A, l, u, xmin, xmax, x0, opt)
%QPS_CLP  Quadratic Program Solver based on CLP - COIN-OR Linear Programming.
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QPS_CLP(H, C, A, L, U, XMIN, XMAX, X0, OPT)
%   A wrapper function providing a MATPOWER standardized interface for using
%   CLP to solve the following QP (quadratic programming) problem:
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
%       X0 : optional starting value of optimization vector X (NOT USED)
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           clp_opt - options struct for CLP, value in
%               verbose overrides these options
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : exit flag, 1 - optimal, -1 - infeasible, -2 - unbounded
%                   -3 - max iterations/time exceeded
%       OUTPUT : struct with fields
%           exitflag - raw CLP exit flag: 0 - optimal, 1 - infeasible,
%               2 - unbounded, 3 - max iterations/time exceeded
%           status - string with explanation of exitflag
%           (iter - depending on build of solver this may contain
%                   the number of iterations)
%       LAMBDA : struct containing the Langrange and Kuhn-Tucker
%           multipliers on the constraints, with fields:
%           mu_l - lower (left-hand) limit on linear constraints
%           mu_u - upper (right-hand) limit on linear constraints
%           lower - lower bound on optimization variables
%           upper - upper bound on optimization variables
%
%   Note the calling syntax is almost identical to that of CLP. The main
%   difference is that the linear constraints are specified with A, L, U
%   instead of A, B, Aeq, Beq.
%
%   Calling syntax options:
%       [x, f, exitflag, output, lambda] = ...
%           qps_clp(H, c, A, l, u, xmin, xmax, x0, opt)
%
%       x = qps_clp(H, c, A, l, u)
%       x = qps_clp(H, c, A, l, u, xmin, xmax)
%       x = qps_clp(H, c, A, l, u, xmin, xmax, x0)
%       x = qps_clp(H, c, A, l, u, xmin, xmax, x0, opt)
%       x = qps_clp(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = qps_clp(...)
%       [x, f] = qps_clp(...)
%       [x, f, exitflag] = qps_clp(...)
%       [x, f, exitflag, output] = qps_clp(...)
%       [x, f, exitflag, output, lambda] = qps_clp(...)
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
%       [x, f, s, out, lambda] = qps_clp(H, c, A, l, u, xmin, [], x0, opt);
%
%   See also CLP.

%   MATPOWER
%   Copyright (c) 2010-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% check for Optimization Toolbox
% if ~have_fcn('quadprog')
%     error('qps_clp: requires the MEX interface to CLP');
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
        error('qps_clp: LP problem must include constraints or variable bounds');
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
if ~issparse(A)
    A = sparse(A);
end
if ~issparse(H)
    H = sparse(H);
end


%% default options
if ~isempty(opt) && isfield(opt, 'verbose') && ~isempty(opt.verbose)
    verbose = opt.verbose;
else
    verbose = 0;
end

%% set up options struct for CLP
if ~isempty(opt) && isfield(opt, 'clp_opt') && ~isempty(opt.clp_opt)
    clp_opt = clp_options(opt.clp_opt);
else
    clp_opt = clp_options;
end

if have_fcn('opti_clp')     %% use OPTI Toolbox verision's MEX interface
    clp_opt.display = verbose;

    [x, f, exitflag, iter, lam] = clp(tril(H), c, A, l, u, xmin, xmax, clp_opt);
    
    output.iter = iter;

    %% repackage lambdas
%     if isempty(x)
%         x = NaN(nx, 1);
%         f = NaN;
%     end
%     if isempty(lam)
%         lambda = struct( ...
%             'mu_l', zeros(nA, 1), ...
%             'mu_u', zeros(nA, 1), ...
%             'lower', zeros(nx, 1), ...
%             'upper', zeros(nx, 1) ...
%         );
%     else
        mu_l = lam.dual_row;
        mu_u = -lam.dual_row;
        lower = lam.dual_col;
        upper = -lam.dual_col;
        mu_l(mu_l < 0) = 0;
        mu_u(mu_u < 0) = 0;
        lower(lower < 0) = 0;
        upper(upper < 0) = 0;
        
        lambda = struct( ...
            'mu_l', mu_l, ...
            'mu_u', mu_u, ...
            'lower', lower, ...
            'upper', upper ...
        );
%     end
else
    clp_opt.verbose = verbose;

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

    %% call the solver
    [x, z, exitflag] = ...
        clp(H, c, Ai, bi, Ae, be, xmin, xmax, clp_opt);

    %% repackage lambdas
    if isempty(x)
        x = NaN(nx, 1);
        f = NaN;
    else
        if isempty(H) || ~any(any(H))
            f = c'*x;
        else
            f = 0.5 * x'*H*x + c'*x;
        end
    end
    if isempty(z)
        lambda = struct( ...
            'mu_l', zeros(nA, 1), ...
            'mu_u', zeros(nA, 1), ...
            'lower', zeros(nx, 1), ...
            'upper', zeros(nx, 1) ...
        );
    else
        neq = length(be);
        nie = length(bi);
        lam.eqlin = z(1:neq);
        lam.ineqlin = z(neq+(1:nie));
    %%-----  MEXCLP DOES NOT RETURN MULTIPLIERS ON VARIABLE BOUNDS :-/  -----
        lam.lower = NaN(nx, 1);
        lam.upper = NaN(nx, 1);
        kl = find(lam.eqlin > 0);   %% lower bound binding
        ku = find(lam.eqlin < 0);   %% upper bound binding

        mu_l = zeros(nA, 1);
        mu_l(ieq(kl)) = lam.eqlin(kl);
        mu_l(igt) = -lam.ineqlin(nlt+(1:ngt));
        mu_l(ibx) = -lam.ineqlin(nlt+ngt+nbx+(1:nbx));

        mu_u = zeros(nA, 1);
        mu_u(ieq(ku)) = -lam.eqlin(ku);
        mu_u(ilt) = -lam.ineqlin(1:nlt);
        mu_u(ibx) = -lam.ineqlin(nlt+ngt+(1:nbx));

        lambda = struct( ...
            'mu_l', mu_l, ...
            'mu_u', mu_u, ...
            'lower', lam.lower, ...
            'upper', lam.upper ...
        );
    end
end

%% set eflag
eflag = -exitflag;
if eflag == 0       %% success
    eflag = 1;
end

%% set status
output.exitflag = exitflag;
switch exitflag
    case 0
        output.status = 'optimal';
    case 1
        output.status = 'primal infeasible';
    case 2
        output.status = 'dual infeasible';
    case 3
        output.status = 'max iterations or time exceeded';
    otherwise
        output.status = 'unknown exit code';
end
