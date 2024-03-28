function [x, f, eflag, output, lambda] = nlps_knitro(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt)
% nlps_knitro - Nonlinear programming (NLP) Solver based on Artelys Knitro.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       NLPS_KNITRO(F_FCN, X0, A, L, U, XMIN, XMAX, GH_FCN, HESS_FCN, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = NLPS_KNITRO(PROBLEM)
%   A wrapper function providing a standardized interface for using
%   Artelys Knitro to solve the following NLP (nonlinear programming) problem:
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
%           knitro_opt  - options struct for Artelys Knitro, value in verbose
%                   overrides these options
%               opts - struct of other values to be passed directly to
%                   Aretelys Knitro via an options file
%               tol_x - termination tol on x
%               tol_f - termination tol on f
%               maxit - maximum number of iterations
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: f_fcn, x0, A, l, u, xmin, xmax,
%                            gh_fcn, hess_fcn, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : exit flag
%           1 = converged
%           negative values = Artelys Knitro specific failure codes
%           (see KNITROMATLAB documentation for details)
%             https://www.artelys.com/docs/knitro/3_referenceManual/knitromatlabReference.html#return-codes-exit-flags
%       OUTPUT : KNITROMATLAB output struct
%           (see KNITROMATLAB documentation for details)
%             https://www.artelys.com/docs/knitro/3_referenceManual/knitromatlabReference.html#output-structure-fields
%           iterations - number of iterations performed
%           funcCount - number of function evaluations
%           constrviolation - maximum of constraint violations
%           firstorderopt - measure of first-order optimality
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
%           nlps_knitro(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
%
%       x = nlps_knitro(f_fcn, x0);
%       x = nlps_knitro(f_fcn, x0, A, l);
%       x = nlps_knitro(f_fcn, x0, A, l, u);
%       x = nlps_knitro(f_fcn, x0, A, l, u, xmin);
%       x = nlps_knitro(f_fcn, x0, A, l, u, xmin, xmax);
%       x = nlps_knitro(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn);
%       x = nlps_knitro(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn);
%       x = nlps_knitro(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
%       x = nlps_knitro(problem);
%               where problem is a struct with fields:
%                   f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt
%                   all fields except 'f_fcn' and 'x0' are optional
%       x = nlps_knitro(...);
%       [x, f] = nlps_knitro(...);
%       [x, f, exitflag] = nlps_knitro(...);
%       [x, f, exitflag, output] = nlps_knitro(...);
%       [x, f, exitflag, output, lambda] = nlps_knitro(...);
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
%       [x, f, exitflag, output, lambda] = nlps_knitro(problem);
%
% See also nlps_master, knitro_nlp, knitromatlab.

%   MP-Opt-Model
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% options
use_ktropts_file = 1;       %% use a Knitro options file to pass options
create_ktropts_file = 0;    %% generate a Knitro options file on the fly

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

%%-----  set up problem  -----
%% build Jacobian and Hessian structure
randx = rand(size(x0));
nonz = 1e-20;
if nonlinear
    [h, g, dhs, dgs] = gh_fcn(randx);
    if issparse(dhs)
        dhs(dhs ~= 0) = nonz;   %% set non-zero entries to tiny value (for adding later)
        dgs(dgs ~= 0) = nonz;   %% set non-zero entries to tiny value (for adding later)
        Js = [dhs'; dgs'];
    else
        dhs = []; dgs = [];
        Js = [];
    end
else
    g = []; h = [];
    dhs = []; dgs = dhs; Js = dhs;
end
neq = length(g);
niq = length(h);
if isempty(hess_fcn)
    [f_, df_, Hs] = f_fcn(randx);   %% cost
else
    lam = struct('eqnonlin', rand(neq, 1), 'ineqnonlin', rand(niq, 1));
    Hs = hess_fcn(randx, lam, 1);
end
if issparse(Hs)
    Hs(Hs ~= 0) = nonz;     %% set non-zero entries to tiny value (for adding later)
else
    Hs = [];
end

%% function definitions
ktr_gh_fcn   = @(x)constraints(x, gh_fcn, dgs, dhs);
ktr_hess_fcn = @(x, lambda)hessian(x, lambda, f_fcn, hess_fcn, Hs);

%% basic optimset options needed for ktrlink
if isfield(opt, 'knitro_opt')
    knitro_opt = opt.knitro_opt;
else
    knitro_opt = struct();
end
kopts = struct();
kopts.GradObj = 'on';
kopts.GradConstr = 'on';
kopts.HessFcn = ktr_hess_fcn;
kopts.Hessian = 'user-supplied';
if issparse(Js)
    kopts.JacobPattern = Js;
end
if issparse(Hs)
    kopts.HessPattern = Hs;
end

if use_ktropts_file
    if isfield(knitro_opt, 'opt_fname') && ~isempty(knitro_opt.opt_fname)
        opt_fname = knitro_opt.opt_fname;
    elseif isfield(knitro_opt, 'opt') && knitro_opt.opt
        opt_fname = sprintf('knitro_user_options_%d.txt', knitro_opt.opt);
    else
        %% create ktropts file
        ktropts.algorithm           = 1;
        ktropts.bar_directinterval  = 0;
        ktropts.outlev              = verbose;
        if isfield(knitro_opt, 'opts')  %% raw Knitro options for options file
            ktropts = nested_struct_copy(ktropts, knitro_opt.opts);
        end
        if isfield(knitro_opt, 'tol_x') && ~isempty(knitro_opt.tol_x)
            ktropts.xtol = knitro_opt.tol_x;
        end
        if isfield(knitro_opt, 'tol_f') && ~isempty(knitro_opt.tol_f)
            ktropts.opttol = knitro_opt.tol_f;
        end
        if isfield(knitro_opt, 'maxit') && knitro_opt.maxit ~= 0
            ktropts.maxit = knitro_opt.maxit;
        end

        opt_fname = write_ktropts(ktropts);
        create_ktropts_file = 1;    %% make a note that I created it
    end
else
    if isfield(knitro_opt, 'opts')  %% raw Knitro options for optimset()
        kopts = nested_struct_copy(kopts, knitro_opt.opts);
    end
    kopts.Algorithm = 'interior-point';
    if isfield(knitro_opt, 'tol_x') && ~isempty(knitro_opt.tol_x)
        kopts.TolX = knitro_opt.tol_x;
    end
    if isfield(knitro_opt, 'tol_f') && ~isempty(knitro_opt.tol_f)
        kopts.TolFun = knitro_opt.tol_f;
    end
    if isfield(knitro_opt, 'maxit') && knitro_opt.maxit ~= 0
        kopts.MaxIter = knitro_opt.maxit;
        kopts.MaxFunEvals = 4 * knitro_opt.maxit;
    end
    if verbose > 1
        kopts.Display = 'iter';
    elseif verbose == 1
        kopts.Display = 'final';
    else
        kopts.Display = 'off';
    end
    opt_fname = [];
end
fmoptions = optimset(kopts);

%%-----  run solver  -----
if have_feature('knitromatlab')
    if have_feature('knitromatlab', 'vnum') < 12.001
        [x, f, eflag, output, Lambda] = knitromatlab(f_fcn, x0, Af, bf, Afeq, bfeq, ...
                                        xmin, xmax, ktr_gh_fcn, [], fmoptions, opt_fname);
    else
        [x, f, eflag, output, Lambda] = knitro_nlp(f_fcn, x0, Af, bf, Afeq, bfeq, ...
                                        xmin, xmax, ktr_gh_fcn, [], fmoptions, opt_fname);
    end
else
    [x, f, eflag, output, Lambda] = ktrlink(f_fcn, x0, Af, bf, Afeq, bfeq, ...
                                    xmin, xmax, ktr_gh_fcn, fmoptions, opt_fname);
end
success = (eflag == 0);
if success
    eflag = 1;      %% success is 1 (not zero), all other values are Knitro return codes
end

%% delete ktropts file
if create_ktropts_file  %% ... but only if I created it
    delete(opt_fname);
end

%% fix Lambdas
%% (shadow prices on equality variable bounds can come back on wrong limit)
kl = find(Lambda.lower > 0 & Lambda.upper == 0);
Lambda.upper(kl) = Lambda.lower(kl);
Lambda.lower(kl) = 0;
ku = find(Lambda.upper < 0 & Lambda.lower == 0);
Lambda.lower(ku) = Lambda.upper(ku);
Lambda.upper(ku) = 0;

%% package up results
[mu_l, mu_u] = convert_lin_constraint_multipliers(Lambda.eqlin, Lambda.ineqlin, ieq, igt, ilt);

lambda = struct( ...
    'lower', -Lambda.lower, ...
    'upper', Lambda.upper, ...
    'eqnonlin', Lambda.eqnonlin, ...
    'ineqnonlin', Lambda.ineqnonlin, ...
    'mu_l', mu_l, ...
    'mu_u', mu_u );


%-----  callback functions  -----
function [h, g, dh, dg] = constraints(x, gh_fcn, dgs, dhs)
if isempty(gh_fcn)
    nx = length(x);
    h = []; g = h; dh = sparse(0, nx); dg = dh;
else
    [h, g, dh, dg] = gh_fcn(x);
    %% add sparse structure (with tiny values) to current matrices to
    %% ensure that sparsity structure matches that supplied
    if issparse(dh)
        dg = dg + dgs;
        dh = dh + dhs;
    end
end

function H = hessian(x, lambda, f_fcn, hess_fcn, Hs)
if isempty(hess_fcn)
    [f_, df_, H] = f_fcn(x);      %% cost
else
    %% add sparse structure (with tiny values) to current matrices to
    %% ensure that sparsity structure matches that supplied
    H = hess_fcn(x, lambda, 1);
    if issparse(H)
        H = H + Hs;
    end
end

%%-----  write_ktropts  -----
function fname = write_ktropts(ktropts)

%% generate file name
fname = sprintf('ktropts_%06d.txt', fix(1e6*rand));

%% open file
[fd, msg] = fopen(fname, 'wt');     %% write options file
if fd == -1
    error('could not create %d : %s', fname, msg);
end

%% write options
fields = fieldnames(ktropts);
for k = 1:length(fields)
    fprintf(fd, '%s %g\n', fields{k}, ktropts.(fields{k}));
end

%% close file
if fd ~= 1
    fclose(fd);
end
