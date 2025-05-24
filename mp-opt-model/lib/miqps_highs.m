function [x, f, eflag, output, lambda] = miqps_highs(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
% miqps_highs - Mixed Integer Linear Program Solver based on HiGHS.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       MIQPS_HIGHS(H, C, A, L, U, XMIN, XMAX, X0, VTYPE, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = MIQPS_HIGHS(PROBLEM)
%   A wrapper function providing a standardized interface for using
%   HiGHS to solve the following MILP (mixed integer linear programming)
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
%           for QP problems, which HiGHS does not handle
%       C : vector of linear cost coefficients
%       A, L, U : define the optional linear constraints. Default
%           values for the elements of L and U are -Inf and Inf,
%           respectively.
%       XMIN, XMAX : optional lower and upper bounds on the
%           X variables, defaults are -Inf and Inf, respectively.
%       X0 : optional starting value of optimization vector X (NOT USED)
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
%           skip_prices (0) - flag that specifies whether or not to
%               skip the price computation stage, in which the problem
%               is re-solved for only the continuous variables, with all
%               others being constrained to their solved values
%           price_stage_warn_tol (1e-7) - tolerance on the objective fcn
%               value and primal variable relative match required to avoid
%               mis-match warning message
%           highs_opt - options struct for HiGHS (see
%               https://ergo-code.github.io/HiGHS/dev/options/definitions/),
%               value in verbose overrides these options
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, vtype, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : exit flag
%           1 = optimal
%           0 = infeasible, unbounded, or otherwise (see
%               output.info.model_status_string for details)
%       OUTPUT : HiGHS info struct
%           (see HiGHS documentation for details)
%       LAMBDA : struct containing the Langrange and Kuhn-Tucker
%           multipliers on the constraints, with fields:
%           mu_l - lower (left-hand) limit on linear constraints
%           mu_u - upper (right-hand) limit on linear constraints
%           lower - lower bound on optimization variables
%           upper - upper bound on optimization variables
%
%   Calling syntax options:
%       [x, f, exitflag, output, lambda] = ...
%           miqps_highs(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
%
%       x = miqps_highs(H, c, A, l, u)
%       x = miqps_highs(H, c, A, l, u, xmin, xmax)
%       x = miqps_highs(H, c, A, l, u, xmin, xmax, x0)
%       x = miqps_highs(H, c, A, l, u, xmin, xmax, x0, vtype)
%       x = miqps_highs(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
%       x = miqps_highs(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, vtype, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = miqps_highs(...)
%       [x, f] = miqps_highs(...)
%       [x, f, exitflag] = miqps_highs(...)
%       [x, f, exitflag, output] = miqps_highs(...)
%       [x, f, exitflag, output, lambda] = miqps_highs(...)
%
%
%   Example: (problem from https://v8doc.sas.com/sashtml/iml/chap8/sect12.htm)
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
%       [x, f, s, out, lambda] = miqps_highs(H, c, A, l, u, xmin, [], x0, vtype, opt);
%
% See also miqps_master, callhighs.

%   MP-Opt-Model
%   Copyright (c) 2010-2025, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% check for HiGHS
% if ~have_feature('highs')
%     error('qps_highs: requires HiGHS (https://highs.dev)');
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
if ~nnz(H)
    isLP = 1;   %% it's an LP
    if isempty(A) && isempty(xmin) && isempty(xmax)
        error('miqps_highs: LP problem must include constraints or variable bounds');
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
    isLP = 0;   %% nope, it's a QP
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
if nA == 0      %% unconstrained
    %% add single non-binding constraint
    A = sparse(1, nx);
    l = -Inf;
    u =  Inf;
end

%% default options
if ~isempty(opt) && isfield(opt, 'verbose') && ~isempty(opt.verbose)
    verbose = opt.verbose;
else
    verbose = 0;
end

%% handle variable types
if isempty(vtype) || isempty(find(vtype == 'B' | vtype == 'I' | ...
        vtype == 'S' | vtype == 'N'))
    mi = 0;
    integrality = [];
else
    mi = 1;
    %% check for (unsupported) MIQP
    if ~isLP
        error('miqps_highs: MIQP problems not supported.');
    end
    %% expand vtype to nx elements if necessary
    if length(vtype) == 1 && nx > 1
        vtype = char(vtype * ones(nx, 1));
    elseif size(vtype, 2) > 1   %% make sure it's a col vector
        vtype = vtype';
    end
    %% clip bounds on 'B' variables to [0, 1]
    k = find(vtype == 'B');
    if ~isempty(k)
        kk = find(xmax(k) > 1);
        xmax(k(kk)) = 1;
        kk = find(xmin(k) < 0);
        xmin(k(kk)) = 0;
    end
    %% form integrality input
    map = struct('C', "c", 'B', "i", 'I', "i", 'S', "sc", 'N', "si");
    integrality = cellfun(@(s) map.(s), cellstr(vtype));
end

%% set up options struct for HiGHS
if ~isempty(opt) && isfield(opt, 'highs_opt') && ~isempty(opt.highs_opt)
    highs_opt = highs_options(opt.highs_opt);
else
    highs_opt = highs_options;
end
if verbose > 1
    highs_opt.output_flag = true;
else
    highs_opt.output_flag = false;
end
highs_opt = highsoptset(highs_opt);

%% get solver type
if isfield(highs_opt, 'solver')
    if ~isLP && ~strcmp(solver, 'choose')
        error('qps_highs: solver = "%s" not supported for QP problems');
%        highs_opt.solver = "choose";
    end
    if mi && ~strcmp(solver, 'choose')
        error('qps_highs: solver = "%s" not supported for MILP problems');
%        highs_opt.solver = "choose";
    end
    solver = highs_opt.solver;
else
    solver = "choose";
end

%% solve model
if verbose
    vn = highsver;
    if isLP
        if mi
            lpqp = 'MILP';
        else
            lpqp = 'LP';
        end
    else
        lpqp = 'QP';
    end
    names = struct( ...
        'choose', 'default', ...
        'simplex', 'dual simplex', ...
        'ipm', 'interior point', ...
        'pdlp', 'primal-dual hybrid gradient' );
    fprintf('HiGHS Version %s -- %s %s solver\n', vn, names.(solver), lpqp);
end
[soln, info, opts, basis] = callhighs(c, A, l, u, xmin, xmax, H, integrality, highs_opt);
if verbose
    fprintf('HiGHS solution status: %s\n', info.model_status_string);
end

%% extract results
x = soln.col_value;
f = info.objective_function_value;
eflag = info.valid && soln.value_valid;
output = info;

%% separate duals into binding lower & upper bounds
if nA == 0      %% unconstrained
    mu_l = [];
    mu_u = [];
else
    mu_l =  soln.row_dual;
    mu_u = -soln.row_dual;
    mu_l(mu_l < 0) = 0;
    mu_u(mu_u < 0) = 0;
end

%% variable bounds
lam_lower =  soln.col_dual;
lam_upper = -soln.col_dual;
lam_lower(lam_lower < 0) = 0;
lam_upper(lam_upper < 0) = 0;

%% package lambdas
lambda = struct( ...
    'mu_l', mu_l, ...
    'mu_u', mu_u, ...
    'lower', lam_lower, ...
    'upper', lam_upper ...
);

if mi && eflag == 1 && (~isfield(opt, 'skip_prices') || ~opt.skip_prices)
    k = find(vtype == 'I' | vtype == 'B' | vtype == 'N' | ...
                (vtype == 'S' & x == 0));
    if length(k) < nx   %% still have some free variables
        if verbose
            fprintf('--- Integer stage complete, starting price computation stage ---\n');
        end
        if isfield(opt, 'price_stage_warn_tol') && ~isempty(opt.price_stage_warn_tol)
            tol = opt.price_stage_warn_tol;
        else
            tol = 1e-7;
        end

        x0 = x;
        x0(k) = round(x0(k));
        xmin(k) = x0(k);
        xmax(k) = x0(k);

        [x_, f_, eflag_, output_, lambda] = qps_highs(H, c, A, l, u, xmin, xmax, x0, opt);
        if eflag ~= eflag_
            error('miqps_highs: EXITFLAG from price computation stage = %d', eflag_);
        end
        if abs(f - f_)/max(abs(f), 1) > tol
            warning('miqps_highs: relative mismatch in objective function value from price computation stage = %g', abs(f - f_)/max(abs(f), 1));
        end
        xn = abs(x);
        xn(xn<1) = 1;
        [mx, k] = max(abs(x - x_) ./ xn);
        if mx > tol
            warning('miqps_highs: max relative mismatch in x from price computation stage = %g (%g)', mx, x(k));
        end
        output.price_stage = output_;
    end
end
