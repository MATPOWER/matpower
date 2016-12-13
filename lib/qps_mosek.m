function [x, f, eflag, output, lambda] = qps_mosek(H, c, A, l, u, xmin, xmax, x0, opt)
%QPS_MOSEK  Quadratic Program Solver based on MOSEK.
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QPS_MOSEK(H, C, A, L, U, XMIN, XMAX, X0, OPT)
%   A wrapper function providing a MATPOWER standardized interface for using
%   MOSEKOPT to solve the following QP (quadratic programming) problem:
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
%           mosek_opt - options struct for MOSEK, value in verbose
%               overrides these options
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : exit flag
%             1 = success
%             0 = terminated at maximum number of iterations
%            -1 = primal or dual infeasible
%           < 0 = the negative of the MOSEK return code
%       OUTPUT : output struct with the following fields:
%           r - MOSEK return code
%           res - MOSEK result struct
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
%           qps_mosek(H, c, A, l, u, xmin, xmax, x0, opt)
%
%       x = qps_mosek(H, c, A, l, u)
%       x = qps_mosek(H, c, A, l, u, xmin, xmax)
%       x = qps_mosek(H, c, A, l, u, xmin, xmax, x0)
%       x = qps_mosek(H, c, A, l, u, xmin, xmax, x0, opt)
%       x = qps_mosek(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = qps_mosek(...)
%       [x, f] = qps_mosek(...)
%       [x, f, exitflag] = qps_mosek(...)
%       [x, f, exitflag, output] = qps_mosek(...)
%       [x, f, exitflag, output, lambda] = qps_mosek(...)
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
%       [x, f, s, out, lambda] = qps_mosek(H, c, A, l, u, xmin, [], x0, opt);
%
%   See also MOSEKOPT.

%   MATPOWER
%   Copyright (c) 2010-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% check for Optimization Toolbox
% if ~have_fcn('mosek')
%     error('qps_mosek: requires MOSEK');
% end

%%----- input argument handling  -----
%% gather inputs
if nargin == 1 && isstruct(H)       %% problem struct
    p = H;
else                                %% individual args
    p = struct('H', H, 'c', c, 'A', A, 'l', l, 'u', u);
    if nargin > 5
        p.xmin = xmin;
        if nargin > 6
            p.xmax = xmax;
            if nargin > 7
                p.x0 = x0;
                if nargin > 8
                    p.opt = opt;
                end
            end
        end
    end
end

%% define nx, set default values for H and c
if ~isfield(p, 'H') || isempty(p.H) || ~any(any(p.H))
    if (~isfield(p, 'A') || isempty(p.A)) && ...
            (~isfield(p, 'xmin') || isempty(p.xmin)) && ...
            (~isfield(p, 'xmax') || isempty(p.xmax))
        error('qps_mosek: LP problem must include constraints or variable bounds');
    else
        if isfield(p, 'A') && ~isempty(p.A)
            nx = size(p.A, 2);
        elseif isfield(p, 'xmin') && ~isempty(p.xmin)
            nx = length(p.xmin);
        else    % if isfield(p, 'xmax') && ~isempty(p.xmax)
            nx = length(p.xmax);
        end
    end
    p.H = sparse(nx, nx);
    qp = 0;
else
    nx = size(p.H, 1);
    qp = 1;
end
if ~isfield(p, 'c') || isempty(p.c)
    p.c = zeros(nx, 1);
end
if ~isfield(p, 'x0') || isempty(p.x0)
    p.x0 = zeros(nx, 1);
end

%% default options
if ~isfield(p, 'opt')
    p.opt = [];
end
if ~isempty(p.opt) && isfield(p.opt, 'verbose') && ~isempty(p.opt.verbose)
    verbose = p.opt.verbose;
else
    verbose = 0;
end
if ~isempty(p.opt) && isfield(p.opt, 'mosek_opt') && ~isempty(p.opt.mosek_opt)
    mosek_opt = mosek_options(p.opt.mosek_opt);
else
    mosek_opt = mosek_options;
end

%% set up problem struct for MOSEK
prob.c = p.c;
if qp
   [prob.qosubi, prob.qosubj, prob.qoval] = find(tril(sparse(p.H)));
end
if isfield(p, 'A') && ~isempty(p.A)
    prob.a = sparse(p.A);
    nA = size(p.A, 1);
else
    nA = 0;
end
if isfield(p, 'l') && ~isempty(p.A)
    prob.blc = p.l;
end
if isfield(p, 'u') && ~isempty(p.A)
    prob.buc = p.u;
end
if isfield(p, 'xmin') && ~isempty(p.xmin)
    prob.blx = p.xmin;
end
if isfield(p, 'xmax') && ~isempty(p.xmax)
    prob.bux = p.xmax;
end

%% A is not allowed to be empty
if ~isfield(prob, 'a') || isempty(prob.a)
    unconstrained = 1;
    prob.a = sparse(1, 1, 1, 1, nx);
    prob.blc = -Inf;
    prob.buc =  Inf;
else
    unconstrained = 0;
end

%%-----  run optimization  -----
if verbose
    s = have_fcn('mosek', 'all');
    if s.vnum < 7
        alg_names = {           %% version 6.x
            'default',              %%  0 : MSK_OPTIMIZER_FREE
            'interior point',       %%  1 : MSK_OPTIMIZER_INTPNT
            '<conic>',              %%  2 : MSK_OPTIMIZER_CONIC
            '<qcone>',              %%  3 : MSK_OPTIMIZER_QCONE
            'primal simplex',       %%  4 : MSK_OPTIMIZER_PRIMAL_SIMPLEX
            'dual simplex',         %%  5 : MSK_OPTIMIZER_DUAL_SIMPLEX
            'primal dual simplex',  %%  6 : MSK_OPTIMIZER_PRIMAL_DUAL_SIMPLEX
            'automatic simplex',    %%  7 : MSK_OPTIMIZER_FREE_SIMPLEX
            '<mixed int>',          %%  8 : MSK_OPTIMIZER_MIXED_INT
            '<nonconvex>',          %%  9 : MSK_OPTIMIZER_NONCONVEX
            'concurrent'            %% 10 : MSK_OPTIMIZER_CONCURRENT
        };
    elseif s.vnum < 8
        alg_names = {           %% version 7.x
            'default',              %%  0 : MSK_OPTIMIZER_FREE
            'interior point',       %%  1 : MSK_OPTIMIZER_INTPNT
            '<conic>',              %%  2 : MSK_OPTIMIZER_CONIC
            'primal simplex',       %%  3 : MSK_OPTIMIZER_PRIMAL_SIMPLEX
            'dual simplex',         %%  4 : MSK_OPTIMIZER_DUAL_SIMPLEX
            'primal dual simplex',  %%  5 : MSK_OPTIMIZER_PRIMAL_DUAL_SIMPLEX
            'automatic simplex',    %%  6 : MSK_OPTIMIZER_FREE_SIMPLEX
            'network simplex',      %%  7 : MSK_OPTIMIZER_NETWORK_PRIMAL_SIMPLEX
            '<mixed int conic>',    %%  8 : MSK_OPTIMIZER_MIXED_INT_CONIC
            '<mixed int>',          %%  9 : MSK_OPTIMIZER_MIXED_INT
            'concurrent',           %% 10 : MSK_OPTIMIZER_CONCURRENT
            '<nonconvex>'           %% 11 : MSK_OPTIMIZER_NONCONVEX
        };
    else
        alg_names = {           %% version 8.x
            '<conic>',              %%  0 : MSK_OPTIMIZER_CONIC
            'dual simplex',         %%  1 : MSK_OPTIMIZER_DUAL_SIMPLEX
            'default',              %%  2 : MSK_OPTIMIZER_FREE
            'automatic simplex',    %%  3 : MSK_OPTIMIZER_FREE_SIMPLEX
            'interior point',       %%  4 : MSK_OPTIMIZER_INTPNT
            '<mixed int>',          %%  5 : MSK_OPTIMIZER_MIXED_INT
            'primal simplex'        %%  6 : MSK_OPTIMIZER_PRIMAL_SIMPLEX
        };
    end
    if qp
        lpqp = 'QP';
    else
        lpqp = 'LP';
    end
    vn = have_fcn('mosek', 'vstr');
    if isempty(vn)
        vn = '<unknown>';
    end
    fprintf('MOSEK Version %s -- %s %s solver\n', ...
            vn, alg_names{mosek_opt.MSK_IPAR_OPTIMIZER+1}, lpqp);
end
cmd = sprintf('minimize echo(%d)', verbose);
[r, res] = mosekopt(cmd, prob, mosek_opt);

%%-----  repackage results  -----
if isfield(res, 'sol')
    if isfield(res.sol, 'bas')
        sol = res.sol.bas;
    else
        sol = res.sol.itr;
    end
    x = sol.xx;
else
    sol = [];
    x = NaN(nx, 1);
end

%%-----  process return codes  -----
if isfield(res, 'symbcon')
    sc = res.symbcon;
else    
    sc = mosek_symbcon;
end
eflag = -r;
msg = '';
switch (r)
    case sc.MSK_RES_OK
        if ~isempty(sol)
%            if sol.solsta == sc.MSK_SOL_STA_OPTIMAL
            if strcmp(sol.solsta, 'OPTIMAL')
                msg = 'The solution is optimal.';
                eflag = 1;
            else
                eflag = -1;
%                 if sol.prosta == sc.MSK_PRO_STA_PRIM_INFEAS
                if strcmp(sol.prosta, 'PRIMAL_INFEASIBLE')
                    msg = 'The problem is primal infeasible.';
%                 elseif sol.prosta == sc.MSK_PRO_STA_DUAL_INFEAS
                elseif strcmp(sol.prosta, 'DUAL_INFEASIBLE')
                    msg = 'The problem is dual infeasible.';
                else
                    msg = sol.solsta;
                end
            end
        end
    case sc.MSK_RES_TRM_MAX_ITERATIONS
        eflag = 0;
        msg = 'The optimizer terminated at the maximum number of iterations.';
    otherwise
        if isfield(res, 'rmsg') && isfield(res, 'rcodestr')
            msg = sprintf('%s : %s', res.rcodestr, res.rmsg);
        else
            msg = sprintf('MOSEK return code = %d', r);
        end
end

if (verbose || r == sc.MSK_RES_ERR_LICENSE || ...
        r == sc.MSK_RES_ERR_LICENSE_EXPIRED || ...
        r == sc.MSK_RES_ERR_LICENSE_VERSION || ...
        r == sc.MSK_RES_ERR_LICENSE_NO_SERVER_SUPPORT || ...
        r == sc.MSK_RES_ERR_LICENSE_FEATURE || ...
        r == sc.MSK_RES_ERR_LICENSE_INVALID_HOSTID || ...
        r == sc.MSK_RES_ERR_LICENSE_SERVER_VERSION || ...
        r == sc.MSK_RES_ERR_MISSING_LICENSE_FILE) ...
        && ~isempty(msg)  %% always alert user of license problems
    fprintf('%s\n', msg);
end

%%-----  repackage results  -----
if nargout > 1
    if r == 0
        f = p.c' * x;
        if ~isempty(p.H)
            f = 0.5 * x' * p.H * x + f;
        end
    else
        f = [];
    end
    if nargout > 3
        output.r = r;
        output.res = res;
        if nargout > 4
            if ~isempty(sol)
                if isfield(sol, 'slx')
                    lambda.lower = sol.slx;
                else
                    lambda.lower = [];
                end
                if isfield(sol, 'sux')
                    lambda.upper = sol.sux;
                else
                    lambda.upper = [];
                end
                if isfield(sol, 'slc')
                    lambda.mu_l  = sol.slc;
                else
                    lambda.mu_l  = [];
                end
                if isfield(sol, 'suc')
                    lambda.mu_u  = sol.suc;
                else
                    lambda.mu_u  = [];
                end
            else
                if isfield(p, 'xmin') && ~isempty(p.xmin)
                    lambda.lower = NaN(nx, 1);
                else
                    lambda.lower = [];
                end
                if isfield(p, 'xmax') && ~isempty(p.xmax)
                    lambda.upper = NaN(nx, 1);
                else
                    lambda.upper = [];
                end
                lambda.mu_l = NaN(nA, 1);
                lambda.mu_u = NaN(nA, 1);
            end
            if unconstrained
                lambda.mu_l  = [];
                lambda.mu_u  = [];
            end
        end
    end
end
