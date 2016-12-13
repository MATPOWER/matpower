function [x, f, eflag, output, lambda] = miqps_ot(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
%MIQPS_OT  Mixed Integer Linear Program Solver based on INTLINPROG.
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       MIQPS_OT(H, C, A, L, U, XMIN, XMAX, X0, OPT)
%   A wrapper function providing a MATPOWER standardized interface for using
%   QUADPROG or LINPROG from the Optimization Toolbox to solve the
%   following QP (quadratic programming) problem:
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
%               allowed values are 'C' (continuous), 'B' (binary) or
%               'I' (integer).
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
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
%           intlinprog_opt - options struct for INTLINPROG, value in
%               verbose overrides these options
%           linprog_opt - options struct for LINPROG, value in
%               verbose overrides these options
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, vtype, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : QUADPROG/LINPROG exit flag
%           (see QUADPROG and LINPROG documentation for details)
%       OUTPUT : QUADPROG/LINPROG output struct
%           (see QUADPROG and LINPROG documentation for details)
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
%           miqps_ot(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
%
%       x = miqps_ot(H, c, A, l, u)
%       x = miqps_ot(H, c, A, l, u, xmin, xmax)
%       x = miqps_ot(H, c, A, l, u, xmin, xmax, x0)
%       x = miqps_ot(H, c, A, l, u, xmin, xmax, x0, vtype)
%       x = miqps_ot(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
%       x = miqps_ot(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, vtype, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = miqps_ot(...)
%       [x, f] = miqps_ot(...)
%       [x, f, exitflag] = miqps_ot(...)
%       [x, f, exitflag, output] = miqps_ot(...)
%       [x, f, exitflag, output, lambda] = miqps_ot(...)
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
%       [x, f, s, out, lambda] = miqps_ot(H, c, A, l, u, xmin, [], x0, vtype, opt);
%
%   See also INTLINPROG, QUADPROG, LINPROG.

%   MATPOWER
%   Copyright (c) 2010-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% check for Optimization Toolbox
% if ~have_fcn('quadprog')
%     error('miqps_ot: requires the Optimization Toolbox');
% end

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
        error('miqps_ot: LP problem must include constraints or variable bounds');
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
if isempty(H) || ~any(any(H))
    isLP = 1;   %% it's an LP
else
    isLP = 0;   %% nope, it's a QP
end

%% default options
if ~isempty(opt) && isfield(opt, 'verbose') && ~isempty(opt.verbose)
    verbose = opt.verbose;
else
    verbose = 0;
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

%% mixed integer?
if isempty(vtype) || isempty(find(vtype == 'B' | vtype == 'I'))
    mi = 0;
    vtype = repmat('C', 1, nx);
else
    mi = 1;
    %% expand vtype to nx elements if necessary
    if length(vtype) == 1 && nx > 1
        vtype = char(vtype * ones(1, nx));
    end
end
%% check bounds for 'B' variables to make sure they are between 0 and 1
k = find(vtype == 'B');
if ~isempty(k)
    kk = find(xmax(k) > 1);
    xmax(k(kk)) = 1;
    kk = find(xmin(k) < 0);
    xmin(k(kk)) = 0;
end
intcon = find(vtype == 'B' | vtype == 'I');

%% set up options
if verbose > 1
    vrb = 'iter';       %% seems to be same as 'final'
elseif verbose == 1
    vrb = 'final';
else
    vrb = 'off';
end

if isLP
    if mi
        ot_opt = optimoptions('intlinprog');
        if ~isempty(opt) && isfield(opt, 'intlinprog_opt') && ~isempty(opt.intlinprog_opt)
            ot_opt = nested_struct_copy(ot_opt, opt.intlinprog_opt);
        end
    else
        ot_opt = optimoptions('linprog');
        if ~isempty(opt) && isfield(opt, 'linprog_opt') && ~isempty(opt.linprog_opt)
            ot_opt = nested_struct_copy(ot_opt, opt.linprog_opt);
        end
    end
else
    if mi
        error('miqps_ot: MIQP problems not supported.');
    end
    ot_opt = optimoptions('quadprog');
    if have_fcn('quadprog_ls')
        ot_opt.Algorithm = 'interior-point-convex';
    else
        ot_opt.LargeScale = 'off';
    end
    if ~isempty(opt) && isfield(opt, 'quadprog_opt') && ~isempty(opt.quadprog_opt)
        ot_opt = nested_struct_copy(ot_opt, opt.quadprog_opt);
    end
end
ot_opt = optimoptions(ot_opt, 'Display', vrb);

%% call the solver
if isLP
    if mi
        [x, f, eflag, output] = ...
            intlinprog(c, intcon, Ai, bi, Ae, be, xmin, xmax, ot_opt);
        lam = [];
    else
        [x, f, eflag, output, lam] = ...
            linprog(c, Ai, bi, Ae, be, xmin, xmax, x0, ot_opt);
    end
else
    [x, f, eflag, output, lam] = ...
        quadprog(H, c, Ai, bi, Ae, be, xmin, xmax, x0, ot_opt);
end

%% repackage lambdas
if isempty(x)
    x = NaN(nx, 1);
end
if isempty(lam) || (isempty(lam.eqlin) && isempty(lam.ineqlin) && ...
                    isempty(lam.lower) && isempty(lam.upper))
    lambda = struct( ...
        'mu_l', NaN(nA, 1), ...
        'mu_u', NaN(nA, 1), ...
        'lower', NaN(nx, 1), ...
        'upper', NaN(nx, 1) ...
    );
else
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
        'lower', lam.lower(1:nx), ...
        'upper', lam.upper(1:nx) ...
    );
end

if mi && eflag == 1 && (~isfield(opt, 'skip_prices') || ~opt.skip_prices)
    if verbose
        fprintf('--- Integer stage complete, starting price computation stage ---\n');
    end
    if isfield(opt, 'price_stage_warn_tol') && ~isempty(opt.price_stage_warn_tol)
        tol = opt.price_stage_warn_tol;
    else
        tol = 1e-7;
    end
    k = intcon;
    x(k) = round(x(k));
    xmin(k) = x(k);
    xmax(k) = x(k);
    x0 = x;
%     opt.linprog_opt.Algorithm = 'dual-simplex';     %% dual-simplex
    
    [x_, f_, eflag_, output_, lambda] = qps_ot(H, c, A, l, u, xmin, xmax, x0, opt);
%     output
%     output.message
    if eflag ~= eflag_
        error('miqps_ot: EXITFLAG from price computation stage = %d', eflag_);
    end
    if abs(f - f_)/max(abs(f), 1) > tol
        warning('miqps_ot: relative mismatch in objective function value from price computation stage = %g', abs(f - f_)/max(abs(f), 1));
    end
    xn = x;
    xn(abs(xn)<1) = 1;
    [mx, k] = max(abs(x - x_) ./ xn);
    if mx > tol
        warning('miqps_ot: max relative mismatch in x from price computation stage = %g (%g)', mx, x(k));
    end
    output.price_stage = output_;
%     output_
%     output_.message
end
