function [x, f, eflag, output, lambda] = qcqps_gurobi(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt)
% qcqps_gurobi - Quadratically Constrained Quadratic Program Solver based on GUROBI.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QCQPS_GUROBI(H, C, Q, B, LQ, UQ, A, L, U, XMIN, XMAX, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = QCQPS_GUROBI(PROBLEM)
%   A wrapper function providing a standardized interface for using
%   GUROBI to solve the following (possibly non-convex) QCQP (quadratically
%   constrained quadratic programming) problem:
%
%       min 1/2 X'*H*X + C'*X
%        X
%
%   subject to
%
%       LQ(i) <= 1/2 X'*Q{i}*X + B(i,:)*X <= UQ(i), i = 1,2,...,nq
%                           (quadratic constraints)
%       L <= A*X <= U       (linear constraints)
%       XMIN <= X <= XMAX   (variable bounds)
%
%   Inputs (all optional except H, C, Q, B, LQ, and UQ):
%       H : matrix (possibly sparse) of quadratic cost coefficients
%       C : vector of linear cost coefficients
%       Q : nq x 1 cell array of sparse quadratic matrices for quadratic constraints
%       B : matrix (possibly sparse) of linear term of quadratic constraints
%       LQ, UQ: define the lower an upper bounds on the quadratic constraints
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
%               3 = even more verbose progress output
%           grb_opt - options struct for GUROBI, value in verbose
%                   overrides these options
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt
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
%           mu_lq - lower (left-hand) limit on quadratic constraints
%           mu_uq - upper (right-hand) limit on quadratic constraints
%           lower - lower bound on optimization variables
%           upper - upper bound on optimization variables
%
%   Calling syntax options:
%       [x, f, exitflag, output, lambda] = ...
%           qcqps_gurobi(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt)
%
%       x = qcqps_gurobi(H, c, Q, B, lq, uq)
%       x = qcqps_gurobi(H, c, Q, B, lq, uq, A, l, u)
%       x = qcqps_gurobi(H, c, Q, B, lq, uq, A, l, u, xmin, xmax)
%       x = qcqps_gurobi(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0)
%       x = qcqps_gurobi(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt)
%       x = qcqps_gurobi(problem), where problem is a struct with fields:
%                       H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'c', are optional, and problem with
%                       linear costs must include constraints
%       x = qcqps_gurobi(...)
%       [x, f] = qcqps_gurobi(...)
%       [x, f, exitflag] = qcqps_gurobi(...)
%       [x, f, exitflag, output] = qcqps_gurobi(...)
%       [x, f, exitflag, output, lambda] = qcqps_gurobi(...)
%
%   Example: (problem from https://docs.gurobi.com/projects/examples/en/current/examples/matlab/qcp.html)
%       H = [];
%       c = [-1;0;0];
%       Q = { sparse([2 0 0; 0 2 0; 0 0 -2]), ...
%             sparse([2 0 0; 0 0 -2; 0 -2 0]) };
%       B = zeros(2,3);
%       lq = [-Inf;-Inf];
%       uq = [0; 0];
%       A = [1 1 1];
%       l = 1;
%       u = 1;
%       xmin = zeros(3,1);
%       xmax = Inf(3,1);
%       x0 = zeros(3,1);
%       opt = struct('verbose', 2);
%       [x, f, s, out, lambda] = ...
%           qcqps_gurobi(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt);
%
% See also qcqps_master, gurobi_options, gurobi.

%   MP-Opt-Model
%   Copyright (c) 2019-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model..
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
    if isfield(p, 'uq'),    uq = p.uq;      else,   uq = [];    end
    if isfield(p, 'lq'),    lq = p.lq;      else,   lq = [];    end
    if isfield(p, 'B'),     B = p.B;        else,   B = [];     end
    if isfield(p, 'Q'),     Q = p.Q;        else,   Q = {};     end
    if isfield(p, 'c'),     c = p.c;        else,   c = [];     end
    if isfield(p, 'H'),     H = p.H;        else,   H = [];     end
else                                %% individual args
    if nargin < 13
        opt = [];
        if nargin < 12
            x0 = [];
            if nargin < 11
                xmax = [];
                if nargin < 10
                    xmin = [];
                    if nargin < 7
                        A = [];
                        l = [];
                        u = [];
                    end
                end
            end
        end
    end
end

%% define nx, nq, nlin, and set default values for missing optional inputs
if ~isempty(Q)
    [nq, ncolQ] = size(Q);
    if ~iscell(Q) || ncolQ ~= 1
        error('qcqps_gurobi: Q must be column vector cell array.')
    end
    sizeQi  = cell2mat(cellfun(@(x)(size(x)), Q, 'UniformOutput', false));
    if any(sizeQi(:) - sizeQi(1))
        error('qcqps_gurobi: All matrices Q{i}, i=1,...,%d must be square of the same size.', nq)
    end
elseif ~isempty(B)
    nq = size(B, 1);
else
    nq = 0;
end

if ~isempty(H)
    [nrowH, ncolH] = size(H);
    if nrowH ~= ncolH
        error('qcqps_gurobi: H must be a square matrix.')
    end
    nx = nrowH;
else
    if ~isempty(c)
        nx = length(c);
    else
        if nq
            nx = size(Q{1},2);
        else
            if ~isempty(B)
                nx = size(B,2);
            else
                if ~isempty(A)
                    nx = size(A,2);
                else
                    error('qcqps_gurobi: inputs arguments H, c, Q, B, and A can not be all empty.')
                end
            end
        end
    end
end

if isempty(H) && isempty(c) && isempty(B) && isempty(A)
    error('qcqps_gurobi: Problem is incomplete: H, c, B and A can not be all empty.')
end
if isempty(H)
    H = sparse(nx,nx);
end
if isempty(c)
    c = sparse(nx,1);
elseif length(c) ~= nx
    error('qcqps_gurobi: Dimension of c (%d) must be iqual to the number of variables (%d).', length(c), nx)
end

if nq
    if isempty(B)
        B = sparse(nq,nx);
    else
        [nrowB, ncolB] = size(B);
        if nrowB ~= nq || ncolB ~= nx
            error('qcqps_gurobi: Dimension of B (%dx%d) should be number of quad constraints times number of variables (%dx%d).', nrowB, ncolB, nq, nx);
        end
    end
else
    Q = {};
    if ~isempty(B)
        [~, ncolB] = size(B);
        if ncolB ~= nx
            error('qcqps_gurobi: The number of columns of matrix B (%d) must be equal to the number of variables (%d).', ncolB, nx);
        end
    else
        B = sparse(nq,nx);
        if ~isempty(lq) ||  ~isempty(uq)
            error('qcqps_gurobi: No quadratic constraints were found. Inputs lq and uq should be empty.')
        end
    end
end

if ~isempty(A)
    [nlin, ncolA] = size(A);
    if ncolA ~= nx
        error('qcqps_gurobi: The number of columns of matrix A (%d) must be equal to the number of variables (%d).', ncolA, nx);
    end
else
    nlin = 0;
    A = sparse(nlin,nx);
end

if isempty(uq)          %% By default, quadratic inequalities are ...
    uq = Inf(nq, 1);        %% ... unbounded above and ...
elseif length(uq) ~= size(B,1)
    error('qcqps_gurobi: Dimension mismatch between uq, Q, and B.')
end
if isempty(lq)
    lq = -Inf(nq, 1);       %% ... unbounded below.
elseif length(lq) ~= size(B,1)
    error('qcqps_gurobi: Dimension mismatch between lq, Q, and B.')
end
if isempty(u)           %% By default, linear inequalities are ...
    u = Inf(nlin, 1);       %% ... unbounded above and ...
elseif length(u) ~= nlin
    error('qcqps_gurobi: Dimension of u (%d) must be iqual to the number of linear constraints (%d).', length(u), nlin')
end
if isempty(l)
    l = -Inf(nlin, 1);      %% ... unbounded below.
elseif length(l) ~= nlin
    error('qcqps_gurobi: Dimension of l (%d) must be iqual to the number of linear constraints (%d).', length(l), nlin')
end
if isempty(xmin)        %% By default, optimization variables are ...
    xmin = -Inf(nx, 1);     %% ... unbounded below and ...
elseif length(xmin) ~= nx
    error('qcqps_gurobi: Dimension of xmin (%d) must be iqual to the number of variables (%d).', length(xmin), nx')
end
if isempty(xmax)
    xmax = Inf(nx, 1);      %% ... unbounded above.
elseif length(xmax) ~= nx
    error('qcqps_gurobi: Dimension of xmax (%d) must be iqual to the number of variables (%d).', length(xmax), nx')
end
if isempty(x0)
    x0 = zeros(nx, 1);
elseif length(x0) ~= nx
    error('qcqps_gurobi: Dimension of x0 (%d) must be iqual to the number of variables (%d).', length(x0), nx')
end

%% default options
if ~isempty(opt) && isfield(opt, 'verbose') && ~isempty(opt.verbose)
    verbose = opt.verbose;
else
    verbose = 0;
end

%% set up options struct for Gurobi
if ~isempty(opt) && isfield(opt, 'grb_opt') && ~isempty(opt.grb_opt)
    g_opt = gurobi_options(opt.grb_opt);
else
    g_opt = gurobi_options;
end
if verbose > 1
    g_opt.LogToConsole = 1;
    g_opt.OutputFlag = 1;
    if verbose > 2
        g_opt.DisplayInterval = 1;
    else
        g_opt.DisplayInterval = 100;
    end
else
    g_opt.LogToConsole = 0;
    g_opt.OutputFlag = 0;
end
g_opt.NonConvex = 2;
g_opt.QCPDual = 1;      %% Turn on the multipliers for quadratic constraints

if ~issparse(A)
    A = sparse(A);
end
if issparse(c)
    c = full(c);
end

%% split up quadratic constraints
if ~isempty(Q)
    [ieq_quad, igt_quad, ilt_quad, Q_quad, B_quad, d_quad] = ...
        convert_quad_constraint(Q, B, lq, uq);
    %% grab some dimensions
    neq_quad = length(ieq_quad);                       %% number of quadratic equalities
    niq_quad = length(ilt_quad) + length(igt_quad);    %% number of quadratic inequalities
else
    Q_quad = {};
    neq_quad = 0; niq_quad = 0;
end

%% split up linear constraints
if isempty(Q) && ~isempty(B)
    [ieq_lin, igt_lin, ilt_lin, A_lin, b_lin] = convert_lin_constraint([A; B], [l; lq], [u; uq]);
else
    [ieq_lin, igt_lin, ilt_lin, A_lin, b_lin] = convert_lin_constraint(A, l, u);
end
%% grab some dimensions
neq_lin = length(ieq_lin);                         %% number of linear equalities
niq_lin = length(ilt_lin) + length(igt_lin);       %% number of linear inequalities

%% set up model
if ~isempty(Q_quad)
    m.quadcon   = cell2struct([ cellfun(@(x)(0.5*x), Q_quad, 'UniformOutput', false),                               ...
                                mat2cell(reshape(B_quad', prod(size(B_quad)) ,[]), nx*ones(neq_quad+niq_quad,1)),   ...
                                num2cell(d_quad,2),                                                                 ...
                                cellstr(char([double('=')*ones(neq_quad,1); double('<')*ones(niq_quad,1)]))], ...
                              {  'Qc'  , ...
                                 'q'   , ...
                                 'rhs' , ...
                                 'sense'  },                                                                  ...
                              2);
end
m.A         = A_lin;
m.rhs       = b_lin;
m.sense     = char([ double('=')*ones(1,neq_lin) double('<')*ones(1,niq_lin) ]);
m.lb        = xmin;
m.ub        = xmax;
m.obj       = c;

%% Call the solver
isemptyQ = cell2mat(cellfun(@(x)(~nnz(x)), Q_quad, 'UniformOutput', false));
if sum(isemptyQ) == (neq_quad + niq_quad)   %% No quadratic terms in quadratic constraints (linear constraints)
    if ~nnz(H)
        lpqcqp = 'LP';
    else
        lpqcqp = 'QP';
        if ~issparse(H)
            H = sparse(H);
        end
        m.Q = 0.5 * H;
    end
else
    lpqcqp = 'QCQP';
    if nnz(H)
        if ~issparse(H)
            H = sparse(H);
        end
        m.Q = 0.5 * H;
    end
end
if verbose
    alg_names = {
        'automatic',
        'primal simplex',
        'dual simplex',
        'interior point',
        'concurrent',
        'deterministic concurrent',
        'deterministic concurrent simplex'
        };
    vn = gurobiver;
    fprintf('Gurobi Version %s -- %s %s solver\n', ...
        vn, alg_names{g_opt.Method+2}, lpqcqp);
end

results = gurobi(m, g_opt);

%% Check for status of the optimization run and prepare output
switch results.status
    case 'LOADED'           %% 1
        eflag = -1;
    case 'OPTIMAL'          %% 2, optimal solution found
        eflag = 1;
    case 'INFEASIBLE'       %% 3
        eflag = -3;
    case 'INF_OR_UNBD'      %% 4
        eflag = -4;
    case 'UNBOUNDED'        %% 5
        eflag = -5;
    case 'CUTOFF'           %% 6
        eflag = -6;
    case 'ITERATION_LIMIT'  %% 7
        eflag = -7;
    case 'NODE_LIMIT'       %% 8
        eflag = -8;
    case 'TIME_LIMIT'       %% 9
        eflag = -9;
    case 'SOLUTION_LIMIT'   %% 10
        eflag = -10;
    case 'INTERRUPTED'      %% 11
        eflag = -11;
    case 'NUMERIC'          %% 12
        eflag = -12;
    case 'SUBOPTIMAL'       %% 13
        eflag = -13;
    case 'INPROGRESS'       %% 14
        eflag = -14;
    case 'USER_OBJ_LIMIT'   %% 15
        eflag = -15;
    case 'WORK_LIMIT'       %% 15
        eflag = -16;
    case 'MEM_LIMIT'        %% 16
        eflag = -17;
    otherwise
        eflag = 0;
end
if nargout > 3
    output = results;
end

%% check for empty results (in case optimization failed)
if ~isfield(results, 'x') || isempty(results.x)
    x = NaN(nx, 1);
    lam.lower   = NaN(nx, 1);
    lam.upper   = NaN(nx, 1);
else
    x = results.x;
    lam.lower   = zeros(nx, 1);
    lam.upper   = zeros(nx, 1);
end
if ~isfield(results, 'objval') || isempty(results.objval)
    f = NaN;
else
    f = results.objval;
end
if ~isfield(results, 'pi') || isempty(results.pi)
    pi  = NaN(length(m.rhs), 1);
else
    pi  = results.pi;
end
if ~isfield(results, 'rc') || isempty(results.rc)
    rc  = NaN(nx, 1);
else
    rc  = results.rc;
end

kl = find(rc > 0);   %% lower bound binding
ku = find(rc < 0);   %% upper bound binding
lam.lower(kl)   =  rc(kl);
lam.upper(ku)   = -rc(ku);

[mu_l, mu_u] = convert_constraint_multipliers( ...
    -pi(1:neq_lin), -pi(neq_lin+(1:niq_lin)), ieq_lin, igt_lin, ilt_lin);

if ~isempty(Q_quad)
    if ~isfield(results, 'qcpi') || isempty(results.qcpi)
        qcpi = NaN(length(m.quadcon), 1);
        % Current version of Gurobi (11.0.3) does not return multipliers for
        % non-convex qcqp. See
        % https://docs.gurobi.com/projects/optimizer/en/current/reference/attributes/constraintquadratic.html#qcpi
    else
        qcpi = results.qcpi;
    end
    [mu_lq, mu_uq] = convert_constraint_multipliers( ...
        -qcpi(1:neq_quad), -qcpi(neq_quad+(1:niq_quad)), ...
        ieq_quad, igt_quad, ilt_quad);

    lambda = struct( ...
    'mu_l'      , mu_l, ...
    'mu_u'      , mu_u, ...
    'mu_lq'     , mu_lq, ...
    'mu_uq'     , mu_uq, ...
    'lower'     , lam.lower, ...
    'upper'     , lam.upper ...
     );
else
    lambda = struct( ...
    'mu_l'      , mu_l, ...
    'mu_u'      , mu_u, ...
    'lower'     , lam.lower, ...
    'upper'     , lam.upper ...
     );
end