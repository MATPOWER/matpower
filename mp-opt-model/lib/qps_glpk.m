function [x, f, eflag, output, lambda] = qps_glpk(H, c, A, l, u, xmin, xmax, x0, opt)
% qps_glpk - Linear Program Solver based on GLPK - GNU Linear Programming Kit.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QPS_GLPK(H, C, A, L, U, XMIN, XMAX, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = QPS_GLPK(PROBLEM)
%   A wrapper function providing a standardized interface for using
%   GLPK to solve the following LP (linear programming) problem:
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
%       H : IGNORED dummy matrix of quadratic cost coefficients
%           for QP problems, which GLPK does not handle
%       C : vector of linear cost coefficients
%       A, L, U : define the optional linear constraints. Default
%           values for the elements of L and U are -Inf and Inf,
%           respectively.
%       XMIN, XMAX : optional lower and upper bounds on the
%           X variables, defaults are -Inf and Inf, respectively.
%       X0 : optional starting value of optimization vector X (NOT USED)
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           glpk_opt - options struct for GLPK, value in verbose
%                   overrides these options
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: H, c, A, l, u, xmin, xmax, x0, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : exit flag, 1 - optimal, <= 0 - infeasible, unbounded or other
%       OUTPUT : output struct with the following fields:
%           errnum - GLPK errnum output arg
%           status - GKPK status output arg
%           runtime - solver run time in seconds
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
%           qps_glpk([], c, A, l, u, xmin, xmax, x0, opt)
%
%       x = qps_glpk([], c, A, l, u)
%       x = qps_glpk([], c, A, l, u, xmin, xmax)
%       x = qps_glpk([], c, A, l, u, xmin, xmax, x0)
%       x = qps_glpk([], c, A, l, u, xmin, xmax, x0, opt)
%       x = qps_glpk(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = qps_glpk(...)
%       [x, f] = qps_glpk(...)
%       [x, f, exitflag] = qps_glpk(...)
%       [x, f, exitflag, output] = qps_glpk(...)
%       [x, f, exitflag, output, lambda] = qps_glpk(...)
%
%   Example: (based on example from 'doc linprog')
%       c = [-5; -4; -6];
%       A = [ 1  -1   1;
%            -3  -2  -4;
%             3   2   0];
%       l = [-Inf; -42; -Inf];
%       u = [20; Inf; 30];
%       xmin = [0; 0; 0];
%       opt = struct('verbose', 2);
%       [x, f, s, out, lambda] = qps_glpk([], c, A, l, u, xmin, [], [], opt);
%
% See also qps_master, glpk.

%   MP-Opt-Model
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% check for Optimization Toolbox
% if ~have_feature('quadprog')
%     error('qps_glpk: requires the MEX interface to GLPK');
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
        error('qps_glpk: LP problem must include constraints or variable bounds');
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
    error('qps_glpk: GLPK handles only LP problems, not QP problems');
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
vtype = repmat('C', nx, 1);

%% set options struct for GLPK
if ~isempty(opt) && isfield(opt, 'glpk_opt') && ~isempty(opt.glpk_opt)
    glpk_opt = glpk_options(opt.glpk_opt);
else
    glpk_opt = glpk_options;
end
glpk_opt.msglev = verbose;

%% call the solver
t0 = tic;
[x, f, errnum, extra] = ...
    glpk(c, AA, bb, xmin, xmax, ctype, vtype, 1, glpk_opt);
output.runtime = toc(t0);

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
