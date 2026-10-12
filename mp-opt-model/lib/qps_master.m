function [x, f, eflag, output, lambda] = qps_master(H, c, A, l, u, xmin, xmax, x0, opt)
% qps_master - Quadratic Program Solver wrapper function.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QPS_MASTER(H, C, A, L, U, XMIN, XMAX, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = QPS_MASTER(PROBLEM)
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
%           alg ('DEFAULT') - determines which solver to use, can be either
%                   a string (new-style) or a numerical alg code (old-style)
%               'DEFAULT' : (or 0) automatic, first available of Gurobi,
%                       CPLEX, MOSEK, Opt Tbx (if MATLAB), HIGHS, GLPK (LPs
%                       only), BPMPD, MIPS
%               'MIPS'    : (or 200) MIPS, MATPOWER Interior Point Solver
%                        pure MATLAB implementation of a primal-dual
%                        interior point method, if mips_opt.step_control = 1
%                        (or alg=250) it uses MIPS-sc, a step controlled
%                        variant of MIPS
%               'BPMPD'   : (or 100) BPMPD_MEX
%               'CLP'     : CLP
%               'CPLEX'   : (or 500) CPLEX
%               'GLPK'    : GLPK, (LP problems only, i.e. empty H matrix)
%               'GUROBI'  : (or 700) Gurobi
%               'HIGHS'   : HiGHS, https://highs.dev
%               'IPOPT'   : (or 400) IPOPT, requires MEX interface to IPOPT
%                           solver, https://github.com/coin-or/Ipopt
%               'MOSEK'   : (or 600) MOSEK
%               'OSQP'    : OSQP, https://osqp.org
%               'OT'      : (or 300) Optimization Toolbox, QUADPROG or LINPROG
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           lazy ([]) - vector of constraint indices (logical or numeric) of
%               lazy constraints, triggers an iterative cutting plane approach,
%               if not empty; set to 'all' to indicate all constraints are lazy
%           lazy_it_lim (Inf) - lazy constraint iteration limit
%           lazy_mode (0) - selects mode for handling lazy constraints
%               0 = add all lazy constraints at iteration lazy_it_lim if any
%                   are still not satisfied
%               1 = scale lazy_thresh by lazy_thresh_multiplier at each
%                   iteration beginning at iteration lazy_it_lim
%               2 = scale lazy_thresh by lazy_thresh_multiplier at each
%                   iteration *beginning when* number of violations increases
%                   from previous iteration, plus lazy_mode = 0 behavior
%               3 = scale lazy_thresh by lazy_thresh_multiplier at each
%                   iteration *in which* number of violations increases
%                   from previous iteration, plus lazy_mode = 0 behavior
%           lazy_thresh ([]) - vector of constraint thresholds for including
%               lazy constraints, size must match total number of constraints
%               or number of lazy constraints
%           lazy_thresh_multiplier (2) - used to scale lazy_thresh as specified
%                by lazy_mode
%           lazy_violation_cost (1) - cost (scalar) of violating a lazy
%               constraint when evaluating violations and thresholds
%           bp_opt      - options vector for BP
%           clp_opt     - options vector for CLP
%           cplex_opt   - options struct for CPLEX
%           glpk_opt    - options struct for GLPK
%           grb_opt     - options struct for GUROBI
%           highs_opt   - options struct for HIGHS
%           ipopt_opt   - options struct for IPOPT
%           knitro_opt  - options struct for KNITRO
%           linprog_opt - options struct for LINPROG
%           mips_opt    - options struct for QPS_MIPS
%           mosek_opt   - options struct for MOSEK
%           osqp_opt    - options struct for OSQP
%           quadprog_opt - options struct for QUADPROG
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : exit flag
%           1 = converged
%           0 or negative values = solver specific failure codes
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
%           qps_master(H, c, A, l, u, xmin, xmax, x0, opt)
%
%       x = qps_master(H, c, A, l, u)
%       x = qps_master(H, c, A, l, u, xmin, xmax)
%       x = qps_master(H, c, A, l, u, xmin, xmax, x0)
%       x = qps_master(H, c, A, l, u, xmin, xmax, x0, opt)
%       x = qps_master(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = qps_master(...)
%       [x, f] = qps_master(...)
%       [x, f, exitflag] = qps_master(...)
%       [x, f, exitflag, output] = qps_master(...)
%       [x, f, exitflag, output, lambda] = qps_master(...)
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
%       [x, f, s, out, lambda] = qps_master(H, c, A, l, u, xmin, [], x0, opt);

%   MP-Opt-Model
%   Copyright (c) 2010-2026, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

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
    %% convert integer codes to string values
    if ~ischar(alg)
        switch alg
            case 0
                alg = 'DEFAULT';
            case 100
                alg = 'BPMPD';
            case 200
                alg = 'MIPS';
                opt.mips_opt.step_control = 0;
            case 250
                alg = 'MIPS';
                opt.mips_opt.step_control = 1;
            case 300
                alg = 'OT';
            case 400
                alg = 'IPOPT';
            case 500
                alg = 'CPLEX';
            case 600
                alg = 'MOSEK';
            case 700
                alg = 'GUROBI';
            case 800
                alg = 'KNITRO';
            otherwise
                error('qps_master: %d is not a valid algorithm code', alg);
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
if ~isempty(opt) && isfield(opt, 'lazy') && any(opt.lazy)
    lazy = opt.lazy;
    [m, n] = size(A);
    %% ensure that lazy is logical
    if ischar(lazy)
        if strcmp(lazy, 'all')
            lazy = true(m, 1);
        else
            error('qps_master: opt.lazy must be a constraint index vector or ''all''');
        end
    elseif ~islogical(lazy)
        if all(lazy >= 1) && all(lazy <= m)
            tmp = false(m, 1);
            tmp(lazy) = true;
            lazy = tmp;
        else
            error('qps_master: if opt.lazy is a numerical vector, all elements must be >= 1 and <= %d', m);
        end
    elseif length(lazy) ~= m
        error('qps_master: if opt.lazy is a logical vector, its length (%d) must match the number of constraints (%d)', length(lazy), m);
    end
    if ~isempty(opt) && isfield(opt, 'lazy_thresh') && ~isempty(opt.lazy_thresh)
        lazy_thresh = opt.lazy_thresh;
        if length(lazy_thresh) ~= m
            if length(lazy_thresh) == sum(lazy)
                tmp = zeros(m, 1);
                tmp(lazy) = lazy_thresh;
                lazy_thresh = tmp;
            else
                error('qps_master: size of opt.lazy_thresh (%d) must match number of lazy (%d) or total (%d) constraints', length(lazy_thresh), sum(lazy), m);
            end
        end
    else
        lazy_thresh = [];
    end
    if ~isempty(opt) && isfield(opt, 'lazy_it_lim') && ~isempty(opt.lazy_it_lim)
        lazy_it_lim = opt.lazy_it_lim;
    else
        lazy_it_lim = Inf;
    end
    if ~isempty(opt) && isfield(opt, 'lazy_mode') && ~isempty(opt.lazy_mode)
        lazy_mode = opt.lazy_mode;
    else
        lazy_mode = 0;
    end
    if ~isempty(opt) && isfield(opt, 'lazy_violation_cost') && ~isempty(opt.lazy_violation_cost)
        lazy_violation_cost = opt.lazy_violation_cost;
    else
        lazy_violation_cost = 1;
    end
    if ~isempty(opt) && isfield(opt, 'lazy_thresh_multiplier') && ~isempty(opt.lazy_thresh_multiplier)
        lazy_thresh_multiplier = opt.lazy_thresh_multiplier;
    else
        lazy_thresh_multiplier = 2;
    end
else
    lazy = [];
end
if strcmp(alg, 'DEFAULT')
    if have_feature('gurobi')       %% use Gurobi by default, if available
        alg = 'GUROBI';
    elseif have_feature('cplex')    %% if not, then CPLEX, if available
        alg = 'CPLEX';
    elseif have_feature('mosek')    %% if not, then MOSEK, if available
        alg = 'MOSEK';
    elseif have_feature('quadprog') && have_feature('matlab')   %% if not, then Opt Tbx, if available in MATLAB
        alg = 'OT';
    elseif have_feature('highs')    %% if not, then HiGHS, if available
        alg = 'HIGHS';
    elseif have_feature('knitro')   %% if not, then Artelys Knitro, if available
        alg = 'KNITRO';
    elseif ~nnz(H) && have_feature('glpk')  %% if not, and
        alg = 'GLPK';               %% prob is LP (not QP), then GLPK, if available
    elseif have_feature('bpmpd')    %% if not, then BPMPD_MEX, if available
        alg = 'BPMPD';
    else                            %% otherwise MIPS
        alg = 'MIPS';
    end
end

%%----- call the appropriate solver  -----
if isempty(lazy)
    switch alg
        case 'BPMPD'                    %% use BPMPD_MEX
            [x, f, eflag, output, lambda] = ...
                qps_bpmpd(H, c, A, l, u, xmin, xmax, x0, opt);

            if eflag == -99
                if verbose
                    mp_printf('         Retrying with QPS_MIPS solver ...\n\n');
                end
                %% save (incorrect) solution from BPMPD
                bpmpd = struct('x', x, 'f', f, 'eflag', eflag, ...
                                'output', output, 'lambda', lambda);
                opt.alg = 'MIPS';
                [x, f, eflag, output, lambda] = ...
                    qps_master(H, c, A, l, u, xmin, xmax, x0, opt);
                output.bpmpd = bpmpd;
            end
        case 'CLP'
            [x, f, eflag, output, lambda] = ...
                qps_clp(H, c, A, l, u, xmin, xmax, x0, opt);
        case 'CPLEX'
            [x, f, eflag, output, lambda] = ...
                qps_cplex(H, c, A, l, u, xmin, xmax, x0, opt);
        case 'GLPK'
            [x, f, eflag, output, lambda] = ...
                qps_glpk(H, c, A, l, u, xmin, xmax, x0, opt);
        case 'GUROBI'
            [x, f, eflag, output, lambda] = ...
                qps_gurobi(H, c, A, l, u, xmin, xmax, x0, opt);
        case 'HIGHS'
            [x, f, eflag, output, lambda] = ...
                qps_highs(H, c, A, l, u, xmin, xmax, x0, opt);
        case 'IPOPT'
            [x, f, eflag, output, lambda] = ...
                qps_ipopt(H, c, A, l, u, xmin, xmax, x0, opt);
        case 'MIPS'
            %% set up options
            if ~isempty(opt) && isfield(opt, 'mips_opt') && ~isempty(opt.mips_opt)
                mips_opt = opt.mips_opt;
            else
                mips_opt = [];
            end
            mips_opt.verbose = verbose;

            %% call solver
            [x, f, eflag, output, lambda] = ...
                qps_mips(H, c, A, l, u, xmin, xmax, x0, mips_opt);
        case 'MOSEK'
            [x, f, eflag, output, lambda] = ...
                qps_mosek(H, c, A, l, u, xmin, xmax, x0, opt);
        case 'OSQP'
            [x, f, eflag, output, lambda] = ...
                qps_osqp(H, c, A, l, u, xmin, xmax, x0, opt);
        case 'OT'                   %% use QUADPROG or LINPROG from Opt Tbx ver 2.x+
            [x, f, eflag, output, lambda] = ...
                qps_ot(H, c, A, l, u, xmin, xmax, x0, opt);
        case 'KNITRO'
            [x, f, eflag, output, lambda] = ...
                qps_knitro(H, c, A, l, u, xmin, xmax, x0, opt);
        otherwise
            fcn = ['qps_' lower(alg)];
            if exist([fcn '.m']) == 2
                [x, f, eflag, output, lambda] = ...
                    feval(fcn, H, c, A, l, u, xmin, xmax, x0, opt);
            else
                error('qps_master: ''%s'' is not a valid algorithm code', alg);
            end
    end
else    %% use cutting plain approach for lazy constraints
    %% set up default inputs
    if isempty(u)           %% linear inequalities are ...
        u = Inf(m, 1);      %% ... unbounded above and ...
    end
    if isempty(l)
        l = -Inf(m, 1);     %% ... unbounded below.
    end
    if isempty(xmax)        %% variables are ...
        xmax = Inf(n, 1);   %% ... unbounded above and ...
    end
    if isempty(xmin)
        xmin = -Inf(n, 1);  %% ... unbounded below.
    end
    if isempty(x0)          %% variables are ...
        x0 = max(xmin, 0) + min(xmax, 0);
    end
    if isempty(c)           %% linear costs are ...
        c = zeros(n, 1);    %% ... zero
    end
    ilazy = 0;
    nlazy = sum(lazy);
    active_rows = ~lazy;
    increase_thresh = false;
    current_thresh_multiplier = 1;
    % make all equality constraints active
    out = {};
    opt_no_lazy = opt;
    opt_no_lazy.lazy = [];
    opt_no_lazy.verbose = max(verbose-1, 0);    %% decrease verbose level
    nn = 1;     %% number of active constraints to add at a time
                %% if problem is unbounded
    infeasible = true;
    x = x0;
    while infeasible
        ilazy = ilazy + 1;      %% iteration counter
        %% solve with active constraints
        t0 = tic;
        [x, f, eflag, output, lambda] = ...
            qps_master(H, c, A(active_rows, :), l(active_rows), u(active_rows), ...
                xmin, xmax, x, opt_no_lazy);
        if verbose
            mp_printf('----- %3d : Solved with %d of %d lazy constraints : eflag = %d (%g sec)\n', ...
                ilazy, sum(lazy & active_rows), nlazy, eflag, toc(t0));
        end
        out{ilazy} = output;
        out{ilazy}.eflag = eflag;

        if ~all(active_rows)
            %% partition constraints and variables into active/inactive
            %% [ la ] <= [ Aaa  0  ] [ xa ] <= [ ua ]
            %% [ li ]    [ Aia Aii ] [ xi ]    [ ui ]
            active_cols = full(any(A(active_rows, :), 1));
            mi = sum(~active_rows);     %% number of inactive rows
            ni = sum(~active_cols);     %% number of inactive cols

            %% check for lazy constraint violations and re-solve, if necessary
            violated = false(m, 1);
            if eflag <= 0 || any(isnan(x))  %% unbounded (or other failure)
                %% add next nn inactive constraints
                inactive = find(~active_rows);
                nn = min(nn, length(inactive));
                violated(inactive(1:nn)) = true;
                nn = nn * 2;
                x = x0;
            else    %% find violations
                t0 = tic;
                Aia_xa = A(~active_rows, active_cols) * x(active_cols);
                if ni == 0  %% no "inactive" vars, evaluate violations directly
                    violated(~active_rows) = ...
                        Aia_xa < l(~active_rows) | ...
                        Aia_xa > u(~active_rows);
                    status = sprintf('direct evaluation (%g sec)', toc(t0));
                else
                    ll = l(~active_rows) - Aia_xa;
                    uu = u(~active_rows) - Aia_xa;
                    AA = [ A(~active_rows, ~active_cols) speye(mi) -speye(mi) ];
                    if nnz(H)
                        HH = [ H(~active_cols, ~active_cols) sparse(ni, 2*mi);
                               sparse(2*mi, ni+2*mi) ];
                        cia = c(~active_cols) + ...
                            ( H(~active_cols, active_cols) + ...
                              H(active_cols, ~active_cols)' ) / 2 * x(active_cols);
                    else
                        HH = [];
                        cia = c(~active_cols);
                    end
                    cc = [ cia; lazy_violation_cost * ones(2*mi, 1) ];
                    xxmin = [ xmin(~active_cols); zeros(2*mi, 1) ];
                    xxmax = [ xmax(~active_cols); Inf(2*mi, 1) ];
                    xx0 = [ x0(~active_cols); zeros(2*mi, 1) ];
                    %% solve LP/QP to find violations, allowing any inactive vars to move
                    [xx, ff, eeflag, ooutput, llambda] = ...
                        qps_master(HH, cc, AA, ll, uu, ...
                            xxmin, xxmax, xx0, opt_no_lazy);
                    if verbose
                        if nnz(HH)
                            status = sprintf('QP : eflag = %d (%g sec)', eeflag, toc(t0));
                        else
                            status = sprintf('LP : eflag = %d (%g sec)', eeflag, toc(t0));
                        end
                    end

                    if eeflag == 1
                        x(~active_cols) = xx(1:ni); %% updated inactive variables
                        violated(~active_rows) = llambda.mu_l + llambda.mu_u > 0;
                        f = c' * x;
                        if nnz(H)
                            f = f + 0.5 * x' * H * x;
                        end
                    else
                        error('qps_master: finding violated lazy constraints failed (eflag = %d).', eeflag);
                    end
                end
            end

            nviolated = sum(violated);
            out{ilazy}.nviolated = nviolated;
            if nviolated        %% update active, re-solve
                if ~isempty(lazy_thresh)
                    if ni == 0  %% no "inactive" vars, evaluate slack directly
                        sl = Aia_xa - l(~active_rows);
                        su = u(~active_rows) - Aia_xa;
                    else
                        AAxx = AA * xx;
                        %% find "slack" columns, i.e. any inactive column with
                        %% zero cost and a single non-zero row in A
                        js = (cia == 0 & sum(A(~active_rows, ~active_cols) ~= 0)' == 1);
                        AAp = AA(:, js);
                        AAn = AAp;
                        AAp(AAp < 0) = 0;   %% positive vals in slack cols
                        AAn(AAn > 0) = 0;   %% negative vals in slack cols
                        %% find slack in actual constraint plus any
                        %% slack provided by these "slack" variables
                        sl = AAxx - ll + ...
                            AAp * (xxmax(js) - xx(js)) + ...
                            AAn * (xxmin(js) - xx(js));
                        su = uu - AAxx + ...
                            AAn * (xx(js) - xxmax(js)) + ...
                            AAp * (xx(js) - xxmin(js));
                    end
                    if lazy_mode ~= 1 && ilazy >= lazy_it_lim
                        violated(~active_rows) = true;  %% include all lazy constraints
                        current_thresh_multiplier = Inf;
                    end
                    switch lazy_mode
                    case 1  %% increase thresh every iteration starting at lazy_it_lim
                        if ilazy >= lazy_it_lim
                            increase_thresh = true;
                        end
                    case 2  %% increase thresh every iteration starting at first violation increase
                        if ilazy > 1 && out{ilazy}.nviolated > out{ilazy-1}.nviolated
                            increase_thresh = true;
                        end
                    case 3  %% increase thresh at every violation increase
                        if ilazy > 1 && out{ilazy}.nviolated > out{ilazy-1}.nviolated
                            increase_thresh = true;
                        else
                            increase_thresh = false;
                        end
                    end
                    if increase_thresh
                        lazy_thresh = lazy_thresh * lazy_thresh_multiplier;
                        current_thresh_multiplier = current_thresh_multiplier * lazy_thresh_multiplier;
                    end
                    out{ilazy}.lazy_thresh_multiplier = current_thresh_multiplier;
                    violated(~active_rows) = violated(~active_rows) | ...
                        sl < lazy_thresh(~active_rows) | ...
                        su < lazy_thresh(~active_rows);
                end

                if verbose
                    if eflag == 1
                        mp_printf('-----       %d violation(s) detected via %s\n', ...
                            nviolated, status);
                    end
                    mp_printf('-----       Adding %d lazy constraint(s)\n', sum(violated));
                end
                active_rows = active_rows | violated;
            else                %% done, no more violations
                if verbose
                    mp_printf('-----       No violations detected via %s\n', status);
                end
                infeasible = false;
            end
        else    %% done, no more inactive constraints left
            infeasible = false;
        end
    end

    %% update lambda
    mu_l = lambda.mu_l;
    mu_u = lambda.mu_u;
    lambda.mu_l = zeros(m, 1);
    lambda.mu_u = zeros(m, 1);
    lambda.mu_l(active_rows) = mu_l;
    lambda.mu_u(active_rows) = mu_u;
    if length(out) > 1  %% keep output of previous iterations
        output.lazy_output = out(1:end-1);
    end
    output.active_constraints = active_rows;
end
if ~isfield(output, 'alg') || isempty(output.alg)
    output.alg = alg;
end
