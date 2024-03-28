function [x, f, eflag, output, lambda] = nlps_fmincon(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt)
% nlps_fmincon - Nonlinear programming (NLP) Solver based on FMINCON.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       NLPS_FMINCON(F_FCN, X0, A, L, U, XMIN, XMAX, GH_FCN, HESS_FCN, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = NLPS_FMINCON(PROBLEM)
%   A wrapper function providing a standardized interface for using
%   FMINCON from the Optimization Toolbox to solve the following NLP
%   (nonlinear programming) problem:
%
%   Minimize a function F(X) beginning from a starting point X0, subject to
%   optional linear and nonlinear constraints and variable bounds.
%
%       min F(X)
%        X
%
%   subject to
%
%       G(X) = 0            (nonlinear equalities)
%       H(X) <= 0           (nonlinear inequalities)
%       L <= A*X <= U       (linear constraints)
%       XMIN <= X <= XMAX   (variable bounds)
%
%   Inputs (all optional except F_FCN and X0):
%       F_FCN : handle to function that evaluates the objective function,
%           its gradients and Hessian for a given value of X. If there
%           are nonlinear constraints, the Hessian information is
%           provided by the HESS_FCN function passed in the 9th argument
%           and is not required here. Calling syntax for this function:
%               [F, DF, D2F] = F_FCN(X)
%       X0 : starting value of optimization vector X
%       A, L, U : define the optional linear constraints. Default
%           values for the elements of L and U are -Inf and Inf,
%           respectively.
%       XMIN, XMAX : optional lower and upper bounds on the
%           X variables, defaults are -Inf and Inf, respectively.
%       GH_FCN : handle to function that evaluates the optional
%           nonlinear constraints and their gradients for a given
%           value of X. Calling syntax for this function is:
%               [H, G, DH, DG] = GH_FCN(X)
%           where the columns of DH and DG are the gradients of the
%           corresponding elements of H and G, i.e. DH and DG are
%           transposes of the Jacobians of H and G, respectively.
%       HESS_FCN : handle to function that computes the Hessian of the
%           Lagrangian for given values of X, lambda and mu, where
%           lambda and mu are the multipliers on the equality and
%           inequality constraints, g and h, respectively. The calling
%           syntax for this function is:
%               LXX = HESS_FCN(X, LAM)
%           where lambda = LAM.eqnonlin and mu = LAM.ineqnonlin.
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           fmincon_opt - options struct for FMINCON, value in verbose
%                   overrides these options
%               opts - struct of other values to be passed directly to
%                   OPTIMSET or OPTIMOPTIONS (overridden by fields below)
%               alg - (4) algorithm codes for FMINCON
%                   1 - active-set (not suitable for large problems)
%                   2 - interior-point, w/default 'bfgs' Hessian approx
%                   3 - interior-point, w/ 'lbfgs' Hessian approx
%                   4 - interior-point, w/exact user-supplied Hessian
%                   5 - interior-point, w/Hessian via finite differences
%                   6 - sqp (not suitable for large problems)
%               tol_x - (1e-4) termination tol on x
%               tol_f - (1e-4) termination tol on f
%               max_it - maximum number of iterations
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: f_fcn, x0, A, l, u, xmin, xmax,
%                            gh_fcn, hess_fcn, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : FMINCON exit flag
%           (see FMINCON documentation for details)
%       OUTPUT : FMINCON output struct
%           (see FMINCON documentation for details)
%       LAMBDA : struct containing the Langrange and Kuhn-Tucker
%           multipliers on the constraints, with fields:
%           eqnonlin - nonlinear equality constraints
%           ineqnonlin - nonlinear inequality constraints
%           mu_l - lower (left-hand) limit on linear constraints
%           mu_u - upper (right-hand) limit on linear constraints
%           lower - lower bound on optimization variables
%           upper - upper bound on optimization variables
%
%   Note the calling syntax is almost identical to that of FMINCON from
%   MathWorks' Optimization Toolbox. The main difference is that the linear
%   constraints are specified with A, L, U instead of A, B, Aeq, Beq. The
%   functions for evaluating the objective function, constraints and Hessian
%   are identical.
%
%   Calling syntax options:
%       [x, f, exitflag, output, lambda] = ...
%           nlps_fmincon(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
%
%       x = nlps_fmincon(f_fcn, x0);
%       x = nlps_fmincon(f_fcn, x0, A, l);
%       x = nlps_fmincon(f_fcn, x0, A, l, u);
%       x = nlps_fmincon(f_fcn, x0, A, l, u, xmin);
%       x = nlps_fmincon(f_fcn, x0, A, l, u, xmin, xmax);
%       x = nlps_fmincon(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn);
%       x = nlps_fmincon(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn);
%       x = nlps_fmincon(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
%       x = nlps_fmincon(problem);
%               where problem is a struct with fields:
%                   f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt
%                   all fields except 'f_fcn' and 'x0' are optional
%       x = nlps_fmincon(...);
%       [x, f] = nlps_fmincon(...);
%       [x, f, exitflag] = nlps_fmincon(...);
%       [x, f, exitflag, output] = nlps_fmincon(...);
%       [x, f, exitflag, output, lambda] = nlps_fmincon(...);
%
%   Example: (problem from https://en.wikipedia.org/wiki/Nonlinear_programming)
%       function [f, df, d2f] = f2(x)
%       f = -x(1)*x(2) - x(2)*x(3);
%       if nargout > 1           %% gradient is required
%           df = -[x(2); x(1)+x(3); x(2)];
%           if nargout > 2       %% Hessian is required
%               d2f = -[0 1 0; 1 0 1; 0 1 0];   %% actually not used since
%           end                                 %% 'hess_fcn' is provided
%       end
%       
%       function [h, g, dh, dg] = gh2(x)
%       h = [ 1 -1 1; 1 1 1] * x.^2 + [-2; -10];
%       dh = 2 * [x(1) x(1); -x(2) x(2); x(3) x(3)];
%       g = []; dg = [];
%       
%       function Lxx = hess2(x, lam, cost_mult)
%       if nargin < 3, cost_mult = 1; end
%       mu = lam.ineqnonlin;
%       Lxx = cost_mult * [0 -1 0; -1 0 -1; 0 -1 0] + ...
%               [2*[1 1]*mu 0 0; 0 2*[-1 1]*mu 0; 0 0 2*[1 1]*mu];
%       
%       problem = struct( ...
%           'f_fcn',    @(x)f2(x), ...
%           'gh_fcn',   @(x)gh2(x), ...
%           'hess_fcn', @(x, lam, cost_mult)hess2(x, lam, cost_mult), ...
%           'x0',       [1; 1; 0], ...
%           'opt',      struct('verbose', 2) ...
%       );
%       [x, f, exitflag, output, lambda] = nlps_fmincon(problem);
%
% See also nlps_master, fmincon.

%   MP-Opt-Model
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%%----- input argument handling  -----
%% gather inputs
if nargin == 1 && isstruct(f_fcn)       %% problem struct
    p = f_fcn;
    f_fcn = p.f_fcn;
    x0 = p.x0;
    nx = size(x0, 1);       %% number of optimization variables
    if isfield(p, 'opt'),       opt = p.opt;            else,   opt = [];       end
    if isfield(p, 'hess_fcn'),  hess_fcn = p.hess_fcn;  else,   hess_fcn = '';  end
    if isfield(p, 'gh_fcn'),    gh_fcn = p.gh_fcn;      else,   gh_fcn = '';    end
    if isfield(p, 'xmax'),      xmax = p.xmax;          else,   xmax = [];      end
    if isfield(p, 'xmin'),      xmin = p.xmin;          else,   xmin = [];      end
    if isfield(p, 'u'),         u = p.u;                else,   u = [];         end
    if isfield(p, 'l'),         l = p.l;                else,   l = [];         end
    if isfield(p, 'A'),         A = p.A;                else,   A=sparse(0,nx); end
else                                    %% individual args
    nx = size(x0, 1);       %% number of optimization variables
    if nargin < 10
        opt = [];
        if nargin < 9
            hess_fcn = '';
            if nargin < 8
                gh_fcn = '';
                if nargin < 7
                    xmax = [];
                    if nargin < 6
                        xmin = [];
                        if nargin < 5
                            u = [];
                            if nargin < 4
                                l = [];
                                A = sparse(0,nx);
                            end
                        end
                    end
                end
            end
        end
    end
end

%% set default argument values if missing
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
if isempty(gh_fcn)
    nonlinear = false;          %% no nonlinear constraints present
else
    nonlinear = true;           %% we have some nonlinear constraints
end

%% default options
if ~isempty(opt) && isfield(opt, 'verbose') && ~isempty(opt.verbose)
    verbose = opt.verbose;
else
    verbose = 0;
end

%% split l <= A*x <= u into less than, equal to, greater than, and
%% doubly-bounded sets
[ieq, igt, ilt, Afeq, bfeq, Af, bf] = convert_lin_constraint(A, l, u);

%% basic optimset options needed for fmincon
if isfield(opt, 'fmincon_opt')
    fmincon_opt = opt.fmincon_opt;
else
    fmincon_opt = struct();
end
if isfield(fmincon_opt, 'opts')
    fopts = fmincon_opt.opts;
else
    fopts = struct();
end
if isfield(fmincon_opt, 'tol_x')
    fopts.TolX = fmincon_opt.tol_x;
end
if isfield(fmincon_opt, 'tol_f')
    fopts.TolFun = fmincon_opt.tol_f;
end
if isfield(fmincon_opt, 'max_it') && fmincon_opt.max_it ~= 0
    fopts.MaxIter = fmincon_opt.max_it;
    fopts.MaxFunEvals = 4 * fmincon_opt.max_it;
end
fopts.GradObj = 'on';
fopts.GradConstr = 'on';
if verbose < 1
    fopts.Display = 'off';
elseif verbose < 2
    fopts.Display = 'final';
elseif verbose < 3
    fopts.Display = 'iter';
else
    fopts.Display = 'testing';
end

%% select algorithm
if have_feature('fmincon_ipm')
    %% set default algorithm, if not specified
    if ~isfield(fmincon_opt, 'alg') || fmincon_opt.alg == 0
        fmincon_opt.alg = 4;  %% interior-point, w/ exact user-supplied Hessian
    end
    if fmincon_opt.alg == 4 && isempty(hess_fcn)
        fmincon_opt.alg = 2;  %% interior-point, w/ default 'bfgs' Hessian approx
    end
    switch fmincon_opt.alg
        case 1      %% active-set (does not use sparse matrices, not suitable for large problems)
            fopts.Algorithm   = 'active-set';
            Af = full(Af);
            Afeq = full(Afeq);
        case 2      %% interior-point, w/ default 'bfgs' Hessian approx
            fopts.Algorithm   = 'interior-point';
        case 3      %% interior-point, w/ 'lbfgs' Hessian approx
            fopts.Algorithm   = 'interior-point';
            fopts.Hessian     = 'lbfgs';
        case 4      %% interior-point, w/ exact user-supplied Hessian
            fopts.Algorithm   = 'interior-point';
            fopts.Hessian     = 'user-supplied';
            fmc_hessian = @(x, lambda)hess_fcn(x, lambda, 1);
            fopts.HessFcn     = fmc_hessian;
        case 5      %% interior-point, w/ finite-diff Hessian
            fopts.Algorithm = 'interior-point';
            fopts.Hessian = 'fin-diff-grads';
            fopts.SubproblemAlgorithm = 'cg';
%             fopts.SubProblem = 'cg';
%             fopts.HessianApproximation = 'finite-difference';
        case 6      %% sqp (does not use sparse matrices, not suitable for large problems)
            fopts.Algorithm = 'sqp';
            Af = full(Af);
            Afeq = full(Afeq);
        otherwise
            error('nlps_fmincon: unknown algorithm specified in ''fmincon.alg'' option');
    end
else
    fopts.LargeScale = 'off';
    Af = full(Af);
    Afeq = full(Afeq);
end
if have_feature('optimoptions')
    fmoptions = optimoptions('fmincon');
    fmoptions = nested_struct_copy(fmoptions, fopts);
else
    fmoptions = optimset(fopts);
end
% fmoptions = optimset(fmoptions, 'DerivativeCheck', 'on', 'FinDiffType', 'central', 'FunValCheck', 'on');
% fmoptions = optimset(fmoptions, 'Diagnostics', 'on');

%%-----  run solver  -----
[x, f, eflag, output, Lambda] = ...
  fmincon(f_fcn, x0, Af, bf, Afeq, bfeq, xmin, xmax, gh_fcn, fmoptions);
success = (eflag > 0);

%% fix Lambdas
%% (shadow prices on equality variable bounds can come back on wrong limit)
kl = find(Lambda.lower < 0 & Lambda.upper == 0);
Lambda.upper(kl) = -Lambda.lower(kl);
Lambda.lower(kl) = 0;
ku = find(Lambda.upper < 0 & Lambda.lower == 0);
Lambda.lower(ku) = -Lambda.upper(ku);
Lambda.upper(ku) = 0;

%% package up results
[mu_l, mu_u] = convert_lin_constraint_multipliers(Lambda.eqlin, Lambda.ineqlin, ieq, igt, ilt);

lambda = struct( ...
    'lower', Lambda.lower, ...
    'upper', Lambda.upper, ...
    'eqnonlin', Lambda.eqnonlin, ...
    'ineqnonlin', Lambda.ineqnonlin, ...
    'mu_l', mu_l, ...
    'mu_u', mu_u );
