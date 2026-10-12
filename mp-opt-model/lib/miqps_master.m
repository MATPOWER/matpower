function [x, f, eflag, output, lambda] = miqps_master(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
% miqps_master - Mixed Integer Quadratic Program Solver wrapper function.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       MIQPS_MASTER(H, C, A, L, U, XMIN, XMAX, X0, VTYPE, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = MIQPS_MASTER(PROBLEM)
%   A common wrapper function for various MILP/MIQP solvers.
%   Solves the following MILP/MIQP (mixed integer linear programming/
%   mixed integer quadratic programming) problem:
%
%       min 1/2 X'*H*X + C'*X
%        X
%
%   subject to
%
%       L <= A*X <= U       (linear constraints)
%       XMIN <= X <= XMAX   (variable bounds)
%       X(i) is integer, for i in I (integer variable constraints)
%       X(b) is binary, for b in B (binary variable constraints)
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
%           alg ('DEFAULT') - determines which solver to use, can be either
%                   a string (new-style) or a numerical alg code (old-style)
%               'DEFAULT' : (or 0) automatic, first available of Gurobi,
%                       CPLEX, MOSEK, Opt Tbx (MILPs only), HIGHS (MILPs only),
%                       GLPK (MILPs only)
%               'CPLEX'   : (or 500) CPLEX
%               'GLPK'    : GLPK, (MILP problems only, i.e. empty H matrix)
%               'GUROBI'  : (or 700) Gurobi
%               'HIGHS'   : HiGHS, https://highs.dev
%               'MOSEK'   : (or 600) MOSEK
%               'OT'      : (or 300) Optimization Toolbox, INTLINPROG
%                           (MILP problems only, i.e. empty H matrix)
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           fix_integer (0) - fix integer variables at value in x0, if true
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
%           lazy_mip_gap (1) - n-element vector of values to use for mip_gap
%               parameter for first n iterations of lazy constraint handling,
%               last element defines final MIP gap tolerance
%           mip_gap ([]) - relative MIP gap tolerance, defaults to solver
%               options or solver default
%           mip_gap_abs ([]) - absolute MIP gap tolerance, defaults to solver
%               options or solver default
%           relax_integer (0) - relax integer constraints, if true
%           skip_prices (0) - flag that specifies whether or not to
%               skip the price computation stage, in which the problem
%               is re-solved for only the continuous variables, with all
%               others being constrained to their solved values
%           price_stage_warn_tol (1e-7) - tolerance on the objective fcn
%               value and primal variable relative match required to avoid
%               mis-match warning message
%           cplex_opt - options struct for CPLEX
%           glpk_opt - options struct for GLPK
%           grb_opt - options struct for GUROBI
%           highs_opt - options struct for HIGHS
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
%           0 or negative values = solver specific failure codes
%       OUTPUT : output struct with the following fields:
%           alg - algorithm code of solver used
%           mip_gap - relative MIP gap for solution, if provided by the solver
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
%           miqps_master(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
%
%       x = miqps_master(H, c, A, l, u)
%       x = miqps_master(H, c, A, l, u, xmin, xmax)
%       x = miqps_master(H, c, A, l, u, xmin, xmax, x0)
%       x = miqps_master(H, c, A, l, u, xmin, xmax, x0, vtype)
%       x = miqps_master(H, c, A, l, u, xmin, xmax, x0, vtype, opt)
%       x = miqps_master(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, vtype, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = miqps_master(...)
%       [x, f] = miqps_master(...)
%       [x, f, exitflag] = miqps_master(...)
%       [x, f, exitflag, output] = miqps_master(...)
%       [x, f, exitflag, output, lambda] = miqps_master(...)
%
%   Example: (problem from MOSEK 6.0 Guided Tour, section  7.13.1
%             https://docs.mosek.com/6.0/toolbox/node009.html)
%       c = [-2; -3];
%       A = sparse([195 273; 4 40]);
%       u = [1365; 140];
%       xmax = [4; Inf];
%       vtype = 'I';
%       opt = struct('verbose', 2);
%       p = struct('c', c, 'A', A, 'u', u, 'xmax', xmax, 'vtype', vtype, 'opt', opt);
%       [x, f, s, out, lam] = miqps_master(p);

%   MP-Opt-Model
%   Copyright (c) 2010-2026, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

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
                error('miqps_master: %d is not a valid algorithm code', alg);
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
            error('miqps_master: opt.lazy must be a constraint index vector or ''all''');
        end
    elseif ~islogical(lazy)
        if all(lazy >= 1) && all(lazy <= m)
            tmp = false(m, 1);
            tmp(lazy) = true;
            lazy = tmp;
        else
            error('miqps_master: if opt.lazy is a numerical vector, all elements must be >= 1 and <= %d', m);
        end
    elseif length(lazy) ~= m
        error('miqps_master: if opt.lazy is a logical vector, its length (%d) must match the number of constraints (%d)', length(lazy), m);
    end
    if ~isempty(opt) && isfield(opt, 'lazy_thresh') && ~isempty(opt.lazy_thresh)
        lazy_thresh = opt.lazy_thresh;
        if length(lazy_thresh) ~= m
            if length(lazy_thresh) == sum(lazy)
                tmp = zeros(m, 1);
                tmp(lazy) = lazy_thresh;
                lazy_thresh = tmp;
            else
                error('miqps_master: size of opt.lazy_thresh (%d) must match number of lazy (%d) or total (%d) constraints', length(lazy_thresh), sum(lazy), m);
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
    else
        if ~nnz(H)                  %% if not, and linear objective
            if have_feature('intlinprog')   %% then Optimization Tbx, if available
                alg = 'OT';
            elseif have_feature('highs')    %% if not, then HiGHS, if available
                alg = 'HIGHS';
            elseif have_feature('glpk')     %% if not, and then GLPK, if available
                alg = 'GLPK';
            else
                if isempty(vtype) || all(vtype == 'C')  %% LP
                    if have_feature('knitro')   %% if not, then Artelys Knitro, if available
                        alg = 'KNITRO';
                    elseif have_feature('bpmpd')%% if not, then BPMPD_MEX, if available
                        alg = 'BPMPD';
                    else                        %% otherwise MIPS
                        alg = 'MIPS';
                    end
                else
                    error('miqps_master: no MILP solvers available - requires CPLEX, Gurobi, MOSEK, INTLINPROG, HIGHS, or GLPK');
                end
            end
        else                        %% if not, and quadratic objective
            if isempty(vtype) || all(vtype == 'C')  %% QP
                if have_feature('quadprog') && have_feature('matlab')   %% if not, then Opt Tbx, if available in MATLAB
                    alg = 'OT';
                elseif have_feature('highs')    %% if not, then HiGHS, if available
                    alg = 'HIGHS';
                elseif have_feature('knitro')   %% if not, then Artelys Knitro, if available
                    alg = 'KNITRO';
                elseif have_feature('bpmpd')    %% if not, then BPMPD_MEX, if available
                    alg = 'BPMPD';
                else                            %% otherwise MIPS
                    alg = 'MIPS';
                end
            else                            %% MIQP
                error('miqps_master: no MIQP solvers available - requires CPLEX, Gurobi, or MOSEK');
            end
        end
    end
end

%% handle relax_integer and fix_integer options
done = false;
fix_integer_presolve = true;
if ~isempty(vtype) && (isfield(opt, 'relax_integer') && opt.relax_integer || ...
        isfield(opt, 'fix_integer') && opt.fix_integer)
    n = max([length(x0) length(xmin) length(xmax) size(A, 2)]);
    if length(vtype) == 1   %% expand if necessary
        vtype = char(vtype * ones(1, n));
    end
    j = (vtype == 'B' | vtype == 'I')';
    if isfield(opt, 'fix_integer') && opt.fix_integer
        x0(j) = round(x0(j));
        if fix_integer_presolve
            x = x0;
            Axj = A(:,j) * x(j);
            cc = c(~j);
            if ~nnz(H)
                HH = [];
            else
                HH = H(~j, ~j);
                cc = cc + 0.5 * (H(~j,j)+H(j,~j)') * x(j);
            end
            [x(~j), f, eflag, output, lambda] = ...
                qps_master(HH, cc,  A(:, ~j), l-Axj, u-Axj, ...
                    xmin(~j), xmax(~j), x(~j), opt);
            f = f + c(j)' * x(j);
            if nnz(H)
                f = f + 0.5 * x(j)' * H(j,j) * x(j);
            end
            mu_lower = zeros(size(x));
            mu_upper = zeros(size(x));
            mu_lower(~j) = lambda.lower;
            mu_upper(~j) = lambda.upper;
            lambda.lower = mu_lower;
            lambda.upper = mu_upper;
            done = true;
        else
            %% fix integer variables
            if ~isempty(j)
                %% expand if necessary
                if length(xmin) == 1
                    xmin = xmin * ones(n, 1);
                elseif isempty(xmin)
                    xmin = -Inf(n, 1);
                end
                if length(xmax) == 1
                    xmax = xmax * ones(n, 1);
                elseif isempty(xmax)
                    xmax = Inf(n, 1);
                end
                xmin(j) = x0(j);
                xmax(j) = x0(j);
            end
        end
    end
    if ~done
        vtype(j) = 'C';
    end
end

%%----- call the appropriate solver  -----
if ~done
    if isempty(lazy)
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
            case 'HIGHS'
                [x, f, eflag, output, lambda] = ...
                    miqps_highs(H, c, A, l, u, xmin, xmax, x0, vtype, opt);
            otherwise
                fcn = ['miqps_' lower(alg)];
                if exist([fcn '.m']) == 2
                    [x, f, eflag, output, lambda] = ...
                        feval(fcn, H, c, A, l, u, xmin, xmax, x0, vtype, opt);
                else
                    error('miqps_master: ''%s'' is not a valid algorithm code', alg);
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
        if isempty(vtype) || length(vtype) ~= n && ...
                isempty(find(vtype == 'B' | vtype == 'I' | ...
                    vtype == 'S' | vtype == 'N', 1))
            vtype = char('C' * ones(1, n));
        elseif length(vtype) == 1 && n > 1
            %% expand vtype to n elements if necessary
            vtype = char(vtype * ones(1, n));
        end
        if verbose
            if isempty(find(vtype == 'B' | vtype == 'I' | ...
                    vtype == 'S' | vtype == 'N', 1))
                mixed_int = false;
            else
                mixed_int = true;
            end
        end
        if isfield(opt, 'lazy_mip_gap') && ~isempty(opt.lazy_mip_gap)
            lazy_mip_gap = opt.lazy_mip_gap;
        else
            lazy_mip_gap = [];
        end
        ivars = (vtype ~= 'C'); %% logical index of integer variables
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
        % opt_no_lazy.skip_prices = 1;
        opt_find_violated = opt_no_lazy;
        opt_find_violated.skip_prices = 0;
        nn = 1;     %% number of active constraints to add at a time
                    %% if problem is unbounded
        infeasible = true;
        x = x0;
        while infeasible
            ilazy = ilazy + 1;      %% iteration counter
            if ~isempty(lazy_mip_gap) && ilazy <= length(lazy_mip_gap)
                opt_no_lazy.mip_gap = lazy_mip_gap(ilazy);
            end
            %% solve with active constraints
            t0 = tic;
            [x, f, eflag, output, lambda] = ...
                miqps_master(H, c, A(active_rows, :), l(active_rows), u(active_rows), ...
                    xmin, xmax, x, vtype, opt_no_lazy);
            if verbose
                if isfield(output, 'mip_gap')
                    mip_gap_out = sprintf(', gap = %g', output.mip_gap);
                else
                    mip_gap_out = '';
                end
                if mixed_int && isfield(opt_no_lazy, 'mip_gap')
                    mip_gap_in = sprintf(', gap tol = %g', opt_no_lazy.mip_gap);
                else
                    mip_gap_in = '';
                end
                mp_printf('----- %3d : Solved with %d of %d lazy constraints%s%s : eflag = %d (%g sec)\n', ...
                    ilazy, sum(lazy & active_rows), nlazy, mip_gap_out, ...
                    mip_gap_in, eflag, toc(t0));
            end
            out{ilazy} = output;
            out{ilazy}.eflag = eflag;

            %% enforce integrality of integer vars, to help with warm starts
            x(ivars) = round(x(ivars));

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
                        vvtype = [ vtype(~active_cols) char('C' * ones(1, 2*mi)) ];
                        %% solve LP/QP to find violations, allowing any inactive vars to move
                        [xx, ff, eeflag, ooutput, llambda] = ...
                            miqps_master(HH, cc, AA, ll, uu, ...
                                xxmin, xxmax, xx0, vvtype, opt_find_violated);
                        if verbose
                            if isempty(find(vvtype == 'B' | vvtype == 'I' | ...
                                    vvtype == 'S' | vvtype == 'N', 1))
                                prefix = '';
                            else
                                prefix = 'MI';
                            end
                            if nnz(HH)
                                status = sprintf('%sQP : eflag = %d (%g sec)', prefix, eeflag, toc(t0));
                            else
                                status = sprintf('%sLP : eflag = %d (%g sec)', prefix, eeflag, toc(t0));
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
                            error('miqps_master: finding violated lazy constraints failed (eflag = %d).', eeflag);
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
                    if ~isempty(lazy_mip_gap) && isfield(output, 'mip_gap') && ...
                            output.mip_gap > lazy_mip_gap(end)
                        %% desired final MIP gap tolerance not met,
                        %% re-solve with desired tolerance
                        if verbose
                            mp_printf('-----       Desired relative MIP tolerance (%g) not met (%g)\n', ...
                                lazy_mip_gap(end), output.mip_gap);
                        end
                        lazy_mip_gap(ilazy+1:end) = lazy_mip_gap(end);
                    else
                        infeasible = false;
                    end
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
end
if ~isfield(output, 'alg') || isempty(output.alg)
    output.alg = alg;
end
