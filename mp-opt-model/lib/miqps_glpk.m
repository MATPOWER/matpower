function [x, f, eflag, output, lambda] = miqps_glpk(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
% miqps_glpk - Mixed Integer Linear Program Solver based on GLPK - GNU Linear Programming Kit.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       MIQPS_GLPK(H, C, A, L, U, XMIN, XMAX, X0, VTYPE, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = MIQPS_GLPK(PROBLEM)
%   A wrapper function providing a standardized interface for using
%   GLKP to solve the following MILP (mixed integer linear programming)
%   problem:
%
%       min C'*X
%        X
%
%   subject to
%
%       L <= A*X <= U       (linear constraints)
%       XMIN <= X <= XMAX   (variable bounds)
%
%   Inputs (all optional except H, C, A and L):
%       H : dummy matrix (possibly sparse) of quadratic cost coefficients
%           for QP problems, which GLPK does not handle
%       C : vector of linear cost coefficients
%       A, L, U : define the optional linear constraints. Default
%           values for the elements of L and U are -Inf and Inf,
%           respectively.
%       XMIN, XMAX : optional lower and upper bounds on the
%           X variables, defaults are -Inf and Inf, respectively.
%       X0 : optional starting value of optimization vector X (NOT USED)
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
%           glpk_opt - options struct for GLPK, value in verbose
%                   overrides these options
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, vtype, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : exit flag, 1 - optimal, <= 0 - infeasible, unbounded or other
%       OUTPUT : struct with errnum and status fields from GLPK output args
%       LAMBDA : struct containing the Langrange and Kuhn-Tucker
%           multipliers on the constraints, with fields:
%           mu_l - lower (left-hand) limit on linear constraints
%           mu_u - upper (right-hand) limit on linear constraints
%           lower - lower bound on optimization variables
%           upper - upper bound on optimization variables
%
%   Note the calling syntax is almost identical to that of GLPK. The main
%   difference is that the linear constraints are specified with A, L, U
%   instead of A, B, Aeq, Beq.
%
%   Calling syntax options:
%       [x, f, exitflag, output, lambda] = ...
%           miqps_glpk(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
%
%       x = miqps_glpk(H, c, A, l, u)
%       x = miqps_glpk(H, c, A, l, u, xmin, xmax)
%       x = miqps_glpk(H, c, A, l, u, xmin, xmax, x0)
%       x = miqps_glpk(H, c, A, l, u, xmin, xmax, x0, vtype)
%       x = miqps_glpk(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
%       x = miqps_glpk(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, vtype, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = miqps_glpk(...)
%       [x, f] = miqps_glpk(...)
%       [x, f, exitflag] = miqps_glpk(...)
%       [x, f, exitflag, output] = miqps_glpk(...)
%       [x, f, exitflag, output, lambda] = miqps_glpk(...)
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
%       [x, f, s, out, lambda] = miqps_glpk(H, c, A, l, u, xmin, [], x0, vtype, opt);
%
% See also miqps_master, glpk.

%   MP-Opt-Model
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% check for Optimization Toolbox
% if ~have_feature('quadprog')
%     error('miqps_glpk: requires the MEX interface to GLPK');
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
        error('miqps_glpk: LP problem must include constraints or variable bounds');
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
    error('miqps_glpk: GLPK handles only LP problems, not QP problems');
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

%% split up linear constraints
[ieq, igt, ilt, AA, bb] = convert_lin_constraint(A, l, u);

%% grab some dimensions
neq = length(ieq);                  %% number of equalities
niq = length(ilt) + length(igt);    %% number of inequalities

ctype = [repmat('S', neq, 1); repmat('U', niq, 1)];

if isempty(vtype) || isempty(find(vtype == 'B' | vtype == 'I'))
    mi = 0;
    vtype = repmat('C', nx, 1);
else
    mi = 1;
    %% expand vtype to nx elements if necessary
    if length(vtype) == 1 && nx > 1
        vtype = char(vtype * ones(nx, 1));
    elseif size(vtype, 2) > 1   %% make sure it's a col vector
        vtype = vtype';
    end
end
%% convert 'B' variables to 'I' and clip bounds to [0, 1]
k = find(vtype == 'B');
if ~isempty(k)
    kk = find(xmax(k) > 1);
    xmax(k(kk)) = 1;
    kk = find(xmin(k) < 0);
    xmin(k(kk)) = 0;
    vtype(k) = 'I';
end

%% set options struct for GLPK
if ~isempty(opt) && isfield(opt, 'glpk_opt') && ~isempty(opt.glpk_opt)
    glpk_opt = glpk_options(opt.glpk_opt);
else
    glpk_opt = glpk_options;
end
glpk_opt.msglev = verbose;

%% call the solver
[x, f, errnum, extra] = ...
    glpk(c, AA, bb, xmin, xmax, ctype, vtype, 1, glpk_opt);

%% set exit flag
if isfield(extra, 'status')             %% status found in extra.status
    output.errnum = errnum;
    output.status = extra.status;
    eflag = -errnum;
    if eflag == 0 && extra.status == 5
        eflag = 1;
    end
else                                    %% status found in errnum
    output.errnum = [];
    output.status = errnum;
    if have_feature('octave')
        if errnum == 180 || errnum == 151 || errnum == 171
            eflag = 1;
        else
            eflag = 0;
        end
    else
        if errnum == 5
            eflag = 1;
        else
            eflag = 0;
        end
    end
end

%% repackage lambdas
if isempty(extra) || ~isfield(extra, 'lambda') || isempty(extra.lambda)
    lambda = struct( ...
        'mu_l', zeros(nA, 1), ...
        'mu_u', zeros(nA, 1), ...
        'lower', zeros(nx, 1), ...
        'upper', zeros(nx, 1) ...
    );
else
    lam.lower = extra.redcosts;
    lam.upper = -extra.redcosts;
    lam.lower(lam.lower < 0) = 0;
    lam.upper(lam.upper < 0) = 0;

    [mu_l, mu_u] = convert_lin_constraint_multipliers( ...
        -extra.lambda(1:neq), -extra.lambda(neq+(1:niq)), ieq, igt, ilt);

    lambda = struct( ...
        'mu_l', mu_l, ...
        'mu_u', mu_u, ...
        'lower', lam.lower(1:nx), ...
        'upper', lam.upper(1:nx) ...
    );
end

if mi && eflag == 1 && (~isfield(opt, 'skip_prices') || ~opt.skip_prices)
    k = find(vtype == 'I' | vtype == 'B');
    if length(k) < nx   %% still have some free variables
        if verbose
            fprintf('--- Integer stage complete, starting price computation stage ---\n');
        end
        if isfield(opt, 'price_stage_warn_tol') && ~isempty(opt.price_stage_warn_tol)
            tol = opt.price_stage_warn_tol;
        else
            tol = 1e-7;
        end
        x(k) = round(x(k));
        xmin(k) = x(k);
        xmax(k) = x(k);
        x0 = x;
        opt.glpk_opt.lpsolver = 1;      %% simplex
        opt.glpk_opt.dual = 0;          %% primal simplex
        if have_feature('octave') && have_feature('octave', 'vnum') >= 3.007
            % opt.glpk_opt.dual = 1;      %% primal simplex
            opt.glpk_opt.dual = 2;      %% dual simplex
        end
    
        [x_, f_, eflag_, output_, lambda] = qps_glpk(H, c, A, l, u, xmin, xmax, x0, opt);
        if eflag ~= eflag_
            error('miqps_glpk: EXITFLAG from price computation stage = %d', eflag_);
        end
        if abs(f - f_)/max(abs(f), 1) > tol
            warning('miqps_glpk: relative mismatch in objective function value from price computation stage = %g', abs(f - f_)/max(abs(f), 1));
        end
        xn = abs(x);
        xn(xn<1) = 1;
        [mx, k] = max(abs(x - x_) ./ xn);
        if mx > tol
            warning('miqps_glpk: max relative mismatch in x from price computation stage = %g (%g)', mx, x(k));
        end
        output.price_stage = output_;
    end
end
