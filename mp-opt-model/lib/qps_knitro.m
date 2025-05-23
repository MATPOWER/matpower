function [x, f, eflag, output, lambda] = qps_knitro(H, c, A, l, u, xmin, xmax, x0, opt)
% qps_knitro - Quadratic Program Solver based on Artelys Knitro.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QPS_KNITRO(H, C, A, L, U, XMIN, XMAX, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = QPS_KNITRO(PROBLEM)
%   A wrapper function providing a standardized interface for using
%   Artelys Knitro to solve the following QP (quadratic programming)
%   problem:
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
%               3 = even more verbose progress output
%           knitro_opt - options struct for Artelys Knitro, value in verbose
%                        overrides these options
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : Artelys Knitro exit flag
%           1 = converged
%           (see Artelys Knitro documentation for details)
%       OUTPUT : Artelys Knitro output struct
%           (see Artelys Knitro documentation for details)
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
%           qps_knitro(H, c, A, l, u, xmin, xmax, x0, opt)
%
%       x = qps_knitro(H, c, A, l, u)
%       x = qps_knitro(H, c, A, l, u, xmin, xmax)
%       x = qps_knitro(H, c, A, l, u, xmin, xmax, x0)
%       x = qps_knitro(H, c, A, l, u, xmin, xmax, x0, opt)
%       x = qps_knitro(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = qps_knitro(...)
%       [x, f] = qps_knitro(...)
%       [x, f, exitflag] = qps_knitro(...)
%       [x, f, exitflag, output] = qps_knitro(...)
%       [x, f, exitflag, output, lambda] = qps_knitro(...)
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
%       [x, f, s, out, lambda] = qps_knitro(H, c, A, l, u, xmin, [], x0, opt);
%
% See also qps_master, knitro_qp, knitro_lp

%   MP-Opt-Model
%   Copyright (c) 2010-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%   and Ray Zimmerman, PSERC Cornell
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

%% define nx, set default values for missing optional inputs
if ~nnz(H)
    if isempty(A) && isempty(xmin) && isempty(xmax)
        error('qps_knitro: LP problem must include constraints or variable bounds');
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

%% set up options struct for Artelys Knitro
if ~isempty(opt) && isfield(opt, 'knitro_opt') && ~isempty(opt.knitro_opt)
    knitro_opt = artelys_knitro_options(opt.knitro_opt);
else
    knitro_opt = knitro_options;
end
if verbose > 1
    knitro_opt.outlev = 3;
    if verbose > 2
        knitro_opt.outlev = 4;
    end
else
     knitro_opt.outlev = 0;
end

if ~issparse(A)
    A = sparse(A);
end
if issparse(c)
    c = full(c);
end

%% split up linear constraints
[ieq, igt, ilt, Ae, be, Ai, bi] = convert_lin_constraint(A, l, u);

%% call the solver
if ~nnz(H)
    lpqp = 'LP';
    [x, f, eflag, output, Lambda] = knitro_lp(c, Ai, bi, Ae, be, xmin, xmax, x0, [], knitro_opt);
else
    lpqp = 'QP';
    if ~issparse(H)
        H = sparse(H);
    end
    [x, f, eflag, output, Lambda] = knitro_qp(H, c, Ai, bi, Ae, be, xmin, xmax, x0, [], knitro_opt);
end

if verbose
    vn = knitrover;
    fprintf('Artelys Knitro Version %s -- %s %s solver\n', ...
        vn, output.algorithm, lpqp);
end

[mu_l, mu_u] = convert_constraint_multipliers(Lambda.eqlin, Lambda.ineqlin, ieq, igt, ilt);

if eflag == 0
    eflag = 1;       %% success is 1 (not zero), all other values are Knitro return codes
end
lambda = struct( ...
    'mu_l', mu_l, ...
    'mu_u', mu_u, ...
    'lower', -1*Lambda.lower, ...
    'upper', Lambda.upper ...
);
