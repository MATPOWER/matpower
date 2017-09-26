function [x, f, eflag, output, lambda] = qps_bpmpd(H, c, A, l, u, xmin, xmax, x0, opt)
%QPS_BPMPD  Quadratic Program Solver based on BPMPD_MEX.
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QPS_BPMPD(H, C, A, L, U, XMIN, XMAX, X0, OPT)
%   A wrapper function providing a MATPOWER standardized interface for using
%   BPMPD_MEX (http://www.pserc.cornell.edu/bpmpd/) to solve the
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
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           bp_opt - options vector for BP, value in verbose
%               overrides these options
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : exit flag,
%             1 = optimal solution
%            -1 = suboptimal solution
%            -2 = infeasible primal
%            -3 = infeasible dual
%           -10 = not enough memory
%           -99 = BPMPD bug: returned infeasible solution
%       OUTPUT : output struct with the following fields:
%           message - exit message
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
%           qps_bpmpd(H, c, A, l, u, xmin, xmax, x0, opt)
%
%       x = qps_bpmpd(H, c, A, l, u)
%       x = qps_bpmpd(H, c, A, l, u, xmin, xmax)
%       x = qps_bpmpd(H, c, A, l, u, xmin, xmax, x0)
%       x = qps_bpmpd(H, c, A, l, u, xmin, xmax, x0, opt)
%       x = qps_bpmpd(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = qps_bpmpd(...)
%       [x, f] = qps_bpmpd(...)
%       [x, f, exitflag] = qps_bpmpd(...)
%       [x, f, exitflag, output] = qps_bpmpd(...)
%       [x, f, exitflag, output, lambda] = qps_bpmpd(...)
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
%       [x, f, s, out, lambda] = qps_bpmpd(H, c, A, l, u, xmin, [], x0, opt);
%
%   See also BPMPD_MEX, http://www.pserc.cornell.edu/bpmpd/.

%   MATPOWER
%   Copyright (c) 2010-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% check for BPMPD_MEX
% if ~have_fcn('bpmpd')
%     error('qps_bpmpd: requires BPMPD_MEX (http://www.pserc.cornell.edu/bpmpd/)');
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
        error('qps_bpmpd: LP problem must include constraints or variable bounds');
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
ieq = find( abs(u-l) <= eps );          %% equality
igt = find( u >=  1e10 & l > -1e10 );   %% greater than, unbounded above
ilt = find( l <= -1e10 & u <  1e10 );   %% less than, unbounded below
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

%% set exit flag
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

%% zero out lambdas smaller than a certain tolerance
y(abs(y) < 1e-9) = 0;
w(abs(w) < 1e-9) = 0;

%% necessary for proper operation of constr.m
if eflag == -2              %% infeasible primal
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

%% Note: BPMPD_MEX has a bug which causes it to return an incorrect
%%       (infeasible) solution for some problems.
%% So we need to double-check for feasibility
if eflag > 0
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
        if verbose > 0
            fprintf('\nWARNING: This version of BPMPD_MEX has a bug which caused it to return\n');
            fprintf(  '         an incorrect (infeasible) solution for this particular problem.\n');
        end
        for err_idx = 1:err_cnt
            fprintf('         %s\n', errs{err_idx});
        end
        if bpmpd_bug_fatal
            error('%s\n%s', ...
                'To avoid this bug in BPMPD_MEX you will need', ...
                'to use a different QP solver for this problem.');
        end
        eflag = -99;
        output.message = [output.message '\nBPMPD bug: returned infeasible solution'];
    end
end
