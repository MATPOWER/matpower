function [x, f, eflag, output, lambda] = qps_matpower(H, c, A, l, u, xmin, xmax, x0, opt)
%QPS_MATPOWER  Quadratic Program Solver for MATPOWER
%   A common wrapper function for various QP solvers. 
%   Solves the following QP (quadratic programming) problem:
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
%       qps_matpower(H, c, A, l, u, xmin, xmax, x0, opt)
%
%   x = qps_matpower(H, c, A, l, u)
%   x = qps_matpower(H, c, A, l, u, xmin, xmax)
%   x = qps_matpower(H, c, A, l, u, xmin, xmax, x0)
%   x = qps_matpower(H, c, A, l, u, xmin, xmax, x0, opt)
%   x = qps_matpower(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'f' and 'x0' are optional
%   x = qps_matpower(...)
%   [x, fval] = qps_matpower(...)
%   [x, fval, exitflag] = qps_matpower(...)
%   [x, fval, exitflag, output] = qps_matpower(...)
%   [x, fval, exitflag, output, lambda] = qps_matpower(...)
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
%           alg (0) - determines which solver to use
%                 0 = use first available solver
%               100 = BPMPD_MEX
%               200 = MIPS, Matlab Interior Point Solver
%                     pure Matlab implementation of a primal-dual
%                     interior point method
%               250 = MIPS-sc, a step controlled variant of MIPS
%               300 = Optimization Toolbox, quadprog() or linprog()
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           max_it (0) - maximum number of iterations allowed
%               0 = use algorithm default
%           bp_opt - options vector for bp()
%           mips_opt - options struct for qps_mips()
%           ot_opt - options struct for quadprog()/linprog()
%       problem : The inputs can alternatively be supplied in a single
%           struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, opt
%
%   Outputs:
%       x : solution vector
%       fval : final objective function value
%       exitflag : exit flag
%           1 = converged
%           0 or negative values = algorithm specific failure codes
%       output : algorithm specific output information
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

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2010 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

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
        error('qps_matpower: LP problem must include constraints or variable bounds');
    else
        if ~isempty(A)
            nx = size(A, 2);
        elseif ~isempty(xmin)
            nx = length(xmin);
        else    % if ~isempty(xmax)
            nx = length(xmax);
        end
    end
%     H = sparse(nx, nx);
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
if ~isempty(opt) && isfield(opt, 'alg') && ~isempty(opt.alg)
    alg = opt.alg;
else
    alg = 0;
end

if ~isempty(opt) && isfield(opt, 'max_it') && ~isempty(opt.max_it)
    max_it = opt.max_it;
else
    max_it = 0;
end

if alg == 0
    if have_fcn('bpmpd')    %% use BPMPD_MEX by default if available
        alg = 100;          %% BPMPD_MEX
    else
        alg = 200;          %% MIPS
    end
end

%%----- call the appropriate solver  -----
if alg == 100                       %% use BPMPD_MEX
    %% make sure args are sparse/full as expected by BPMPD
    if ~isempty(H)
        if ~issparse(H)
            H = sparse(H);
        end
    end
    if ~issparse(A)
        A = sparse(A);
    end
    if issparse(c)
        c = full(c);                %% avoid a crash
    end

    %% split up linear constraints
    ieq = find( abs(u-l) <= eps );            %% equality
    igt = find( u >=  1e10 & l > -1e10 );     %% greater than, unbounded above
    ilt = find( l <= -1e10 & u <  1e10 );     %% less than, unbounded below
    ibx = find( (abs(u-l) > eps) & (u < 1e10) & (l > -1e10) );
    Ae = A(ieq, :);
    be = u(ieq);
    Ai  = [ A(ilt, :); -A(igt, :); A(ibx, :); -A(ibx, :) ];
    bi  = [ u(ilt);    -l(igt);    u(ibx);    -l(ibx)];

    %% grab some dimensions
    neq = length(ieq);      %% number of equality constraints
    niq = length(bi);       %% number of inequality constraints
    nlt = length(ilt);      %% number of upper bounded linear inequalities
    ngt = length(igt);      %% number of lower bounded linear inequalities
    nbx = length(ibx);      %% number of doubly bounded linear inequalities

    %% set up linear constraints
    if neq || niq
        AA = [Ae; Ai];
        bb = [be; bi];
        ee = [zeros(neq, 1); -ones(niq, 1)];
    else
        AA = []; bb = []; ee = [];
    end

    %% set up variable bounds and initial value
    if ~isempty(xmin)
        llist = find(xmin > -1e15);  % interpret limits <= -1e15 as unbounded
        if isempty(llist)
            llist = [];
            lval  = [];
        else
            lval = xmin(llist);
        end
    else
        llist = [];
        lval = [];
    end
    if ~isempty(xmax)
        ulist = find(xmax < 1e15);   % interpret limits >= 1e15 as unbounded
        if isempty(ulist)
            ulist = [];
            uval  = [];
        else
            uval = xmax(ulist);
        end
    else
        ulist = [];
        uval = [];
    end

    %% set up options
    if ~isempty(opt) && isfield(opt, 'bp_opt') && ~isempty(opt.bp_opt)
        bp_opt = opt.bp_opt;
    else
        bp_opt = bpopt;         %% use default options
        % bp_opt(14)= 1e-3;   % TPIV1  first relative pivot tolerance (desired)
        % bp_opt(20)= 1e-8;   % TOPT1  stop if feasible and rel. dual gap less than this
        % bp_opt(22)= 1e-7;   % TFEAS1 relative primal feasibility tolerance
        % bp_opt(23)= 1e-7;   % TFEAS2 relative dual feasibility tolerance
        % bp_opt(29)= 1e-9;   % TRESX  acceptable primal residual
        % bp_opt(30)= 1e-9;   % TRESY  acceptable dual residual
        % bp_opt(38)= 2;      % SMETHOD1 prescaling method
    end
    if verbose > 1
        prnlev = 1;
    else
        prnlev = 0;
    end
    if strcmp(computer, 'PCWIN')
        if prnlev
            fprintf('Windows version of BPMPD_MEX cannot print to screen.\n');
        end
        prnlev = 0;   % The Windows incarnation of bp was born mute and deaf,
    end               % probably because of acute shock after realizing its fate.
                      % Can't be allowed to try to speak or its universe crumbles.
    if max_it
        bp_opt(26) = max_it;    %% MAXITER
    end

    %% call the solver
    [x, y, s, w, output.message] = bp(H, AA, bb, c, ee, llist, lval, ...
                                ulist, uval, bp_opt, prnlev);

    %% compute final objective
    if nargout > 1
        f = 0;
        if ~isempty(c)
            f = f + c' * x;
        end
        if ~isempty(H)
            f = f + 0.5 * x' * H * x;
        end
    end

    %% zero out lambdas smaller than a certain tolerance
    y(abs(y) < 1e-9) = 0;
    w(abs(w) < 1e-9) = 0;

    %% necessary for proper operation of constr.m
    if strcmp(output.message, 'infeasible primal')
        y = zeros(size(y));
        w = zeros(size(w));
    end

    %% repackage lambdas
    lam.eqlin   = -y(1:neq);
    lam.ineqlin = -y(neq+(1:niq));
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

    lam.lower = zeros(nx, 1);
    lam.upper = zeros(nx, 1);
    kl = find(w > 0);       %% lower bound binding
    ku = find(w < 0);       %% upper bound binding
    lam.lower(kl) = w(kl);
    lam.upper(ku) = -w(ku);

    lambda = struct( ...
        'mu_l', mu_l, ...
        'mu_u', mu_u, ...
        'lower', lam.lower, ...
        'upper', lam.upper ...
    );

    if strcmp(output.message, 'optimal solution')
        eflag = 1;
    elseif strcmp(output.message, 'suboptimal solution')
        eflag = -1;
    elseif strcmp(output.message, 'infeasible primal')
        eflag = -2;
    elseif strcmp(output.message, 'infeasible dual')
        eflag = -3;
    elseif strcmp(output.message, 'not enough memory')
        eflag = -10;
    else
        eflag = 0;
    end
    if eflag > 0
        %% double-check feasibility
        %% Note: BPMPD_MEX has a bug which causes it to return an incorrect
        %%       (infeasible) solution for some problems.
        bpmpd_bug_fatal = 0;
        err_tol = 5e-4;
        if ~isempty(xmin)
            lb_violation = xmin - x;
            if verbose > 1
                fprintf('max variable lower bound violatation: %g\n', max(lb_violation));
            end
        else
            lb_violation = zeros(nx, 1);
        end
        if ~isempty(xmax)
            ub_violation = x - xmax;
            if verbose > 1
                fprintf('max variable upper bound violation: %g\n', max(ub_violation));
            end
        else
            ub_violation = zeros(nx, 1);
        end
        if neq > 0
            eq_violation = abs( Ae * x - be );
            if verbose > 1
                fprintf('max equality constraint violation: %g\n', max(eq_violation));
            end
        else
            eq_violation = zeros(neq, 1);
        end
        if niq
            ineq_violation = Ai * x - bi;
            if verbose > 1
                fprintf('max inequality constraint violation: %g\n', max(ineq_violation));
            end
        else
            ineq_violation = zeros(niq, 1);
        end
        if any( [ lb_violation;
                  ub_violation;
                  eq_violation;
                  ineq_violation ] > err_tol)
            err_cnt = 0;
            if any( lb_violation > err_tol )
                err_cnt = err_cnt + 1;
                errs{err_cnt} = ...
                    sprintf('variable lower bound violated by %g', ...
                        max(lb_violation));
            end
            if any( ub_violation > err_tol )
                err_cnt = err_cnt + 1;
                errs{err_cnt} = ... 
                    sprintf('variable upper bound violated by %g', ...
                        max(ub_violation));
            end
            if any( eq_violation > err_tol )
                err_cnt = err_cnt + 1;
                errs{err_cnt} = ... 
                    sprintf('equality constraint violated by %g', ...
                        max(eq_violation));
            end
            if any( ineq_violation > err_tol )
                err_cnt = err_cnt + 1;
                errs{err_cnt} = ... 
                    sprintf('inequality constraint violated by %g', ...
                        max(ineq_violation));
            end
            fprintf('\nWARNING: This version of BPMPD_MEX has a bug which caused it to return\n');
            fprintf(  '         an incorrect (infeasible) solution for this particular problem.\n');
            for err_idx = 1:err_cnt
                fprintf('         %s\n', errs{err_idx});
            end
            if bpmpd_bug_fatal
                error('%s\n%s', ...
                    'To avoid this bug in BPMPD_MEX you will need', ...
                    'to use a different QP solver for this problem.');
            end
            eflag = -99;
            output.message = 'BPMPD bug: returned infeasible solution';
        end
    end

    if eflag == -99
        fprintf('         Retrying with MIPS QP solver ...\n\n');
        if ~isempty(opt) && isfield(opt, 'mips_opt') && ~isempty(opt.mips_opt)
            mips_opt = opt.mips_opt;
        else
            mips_opt = [];
        end
        if max_it
            mips_opt.max_it = max_it;
        end
        mips_opt.verbose = verbose;
        [x, f, eflag, output, lambda] = ...
            qps_mips(H, c, A, l, u, xmin, xmax, x0, mips_opt);
    end
elseif alg == 200 || alg == 250     %% use MIPS or sc-MIPS
    %% set up options
    if ~isempty(opt) && isfield(opt, 'mips_opt') && ~isempty(opt.mips_opt)
        mips_opt = opt.mips_opt;
    else
        mips_opt = [];
    end
    if alg == 200
        mips_opt.step_control = 0;
    else
        mips_opt.step_control = 1;
    end
    if max_it
        mips_opt.max_it = max_it;
    end
    mips_opt.verbose = verbose;
    
    mlver = ver('matlab');
    if str2double(mlver.Version(1)) < 7    %% anonymous functions not available
        solver = 'qps_mips6';
    else
        solver = 'qps_mips';
    end

    %% call solver
    [x, f, eflag, output, lambda] = ...
        feval(solver, H, c, A, l, u, xmin, xmax, x0, mips_opt);
elseif alg == 300                   %% use quadprog() or linprog() from Opt Tbx ver 2.x+
    if have_fcn('quadprog')
        %% split up linear constraints
        ieq = find( abs(u-l) <= eps );            %% equality
        igt = find( u >=  1e10 & l > -1e10 );     %% greater than, unbounded above
        ilt = find( l <= -1e10 & u <  1e10 );     %% less than, unbounded below
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
            ot_opt = [];
        end
        if isempty(ot_opt)
            if isempty(H)
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

        %% call solver
        if isempty(H)
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
    else
        error('qps_matpower: alg code %d requires the Optimization Toolbox.', alg);
    end
else
    error('qps_matpower: %d is not a valid algorithm code', alg);
end
