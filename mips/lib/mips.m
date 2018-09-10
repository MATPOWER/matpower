function [x, f, eflag, output, lambda] = mips(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt)
%MIPS  MATPOWER Interior Point Solver.
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       MIPS(F_FCN, X0, A, L, U, XMIN, XMAX, GH_FCN, HESS_FCN, OPT)
%   Primal-dual interior point method for NLP (nonlinear programming).
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
%           linsolver ('') - linear system solver for solving update steps
%               ''          = default solver (currently same as '\')
%               '\'         = built-in \ operator
%               'PARDISO'   = PARDISO solver (if available)
%           feastol (1e-6) - termination tolerance for feasibility
%               condition
%           gradtol (1e-6) - termination tolerance for gradient
%               condition
%           comptol (1e-6) - termination tolerance for complementarity
%               condition
%           costtol (1e-6) - termination tolerance for cost condition
%           max_it (150) - maximum number of iterations
%           step_control (0) - set to 1 to enable step-size control
%           sc.red_it (20) - maximum number of step-size reductions if
%               step-control is on
%           cost_mult (1) - cost multiplier used to scale the objective
%               function for improved conditioning. Note: This value is
%               also passed as the 3rd argument to the Hessian evaluation
%               function so that it can appropriately scale the
%               objective function term in the Hessian of the
%               Lagrangian.
%           xi (0.99995) - constant used in alpha updates*
%           sigma (0.1) - centering parameter*
%           z0 (1) - used to initialize slack variables*
%           alpha_min (1e-8) - returns EXITFLAG = -1 if either alpha
%               parameter becomes smaller than this value*
%           rho_min (0.95) - lower bound on rho_t*
%           rho_max (1.05) - upper bound on rho_t*
%           mu_threshold (1e-5) - KT multipliers smaller than this value
%               for non-binding constraints are forced to zero
%           max_stepsize (1e10) - returns EXITFLAG = -1 if the 2-norm of
%               the reduced Newton step exceeds this value*
%               * see the corresponding Appendix in the manual for details
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: f_fcn, x0, A, l, u, xmin, xmax,
%                            gh_fcn, hess_fcn, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : exit flag
%           1 = first order optimality conditions satisfied
%           0 = maximum number of iterations reached
%           -1 = numerically failed
%       OUTPUT : output struct with fields:
%           iterations - number of iterations performed
%           hist - struct array with trajectories of the following:
%                   feascond, gradcond, compcond, costcond, gamma,
%                   stepsize, obj, alphap, alphad
%           message - exit message
%       LAMBDA : struct containing the Langrange and Kuhn-Tucker
%           multipliers on the constraints, with fields:
%           eqnonlin - nonlinear equality constraints
%           ineqnonlin - nonlinear inequality constraints
%           mu_l - lower (left-hand) limit on linear constraints
%           mu_u - upper (right-hand) limit on linear constraints
%           lower - lower bound on optimization variables
%           upper - upper bound on optimization variables
%
%   Note the calling syntax is almost identical to that of FMINCON
%   from MathWorks' Optimization Toolbox. The main difference is that
%   the linear constraints are specified with A, L, U instead of
%   A, B, Aeq, Beq. The functions for evaluating the objective
%   function, constraints and Hessian are identical.
%
%   Calling syntax options:
%       [x, f, exitflag, output, lambda] = ...
%           mips(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
%
%       x = mips(f_fcn, x0);
%       x = mips(f_fcn, x0, A, l);
%       x = mips(f_fcn, x0, A, l, u);
%       x = mips(f_fcn, x0, A, l, u, xmin);
%       x = mips(f_fcn, x0, A, l, u, xmin, xmax);
%       x = mips(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn);
%       x = mips(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn);
%       x = mips(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
%       x = mips(problem);
%               where problem is a struct with fields:
%                   f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt
%                   all fields except 'f_fcn' and 'x0' are optional
%       x = mips(...);
%       [x, f] = mips(...);
%       [x, f, exitflag] = mips(...);
%       [x, f, exitflag, output] = mips(...);
%       [x, f, exitflag, output, lambda] = mips(...);
%
%   Example: (problem from http://en.wikipedia.org/wiki/Nonlinear_programming)
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
%       [x, f, exitflag, output, lambda] = mips(problem);
%
%   Ported by Ray Zimmerman from C code written by H. Wang for his
%   PhD dissertation:
%     "On the Computation and Application of Multi-period
%     Security-Constrained Optimal Power Flow for Real-time
%     Electricity Market Operations", Cornell University, May 2007.
%
%   See also:
%     H. Wang, C. E. Murillo-Sanchez, R. D. Zimmerman, R. J. Thomas,
%     "On Computational Issues of Market-Based Optimal Power Flow",
%     IEEE Transactions on Power Systems, Vol. 22, No. 3, Aug. 2007,
%     pp. 1185-1193.

%   MIPS
%   Copyright (c) 2009-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MIPS.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mips for more info.

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
    gn = []; hn = [];
else
    nonlinear = true;           %% we have some nonlinear constraints
end

%% default options
if isempty(opt)
    opt = struct('verbose', 0);
end
if ~isfield(opt, 'linsolver') || isempty(opt.linsolver)
    opt.linsolver = '';
end
if ~isfield(opt, 'feastol') || isempty(opt.feastol) || opt.feastol == 0
    opt.feastol = 1e-6;
end
if ~isfield(opt, 'gradtol') || isempty(opt.gradtol) || opt.gradtol == 0
    opt.gradtol = 1e-6;
end
if ~isfield(opt, 'comptol') || isempty(opt.comptol) || opt.comptol == 0
    opt.comptol = 1e-6;
end
if ~isfield(opt, 'costtol') || isempty(opt.costtol) || opt.costtol == 0
    opt.costtol = 1e-6;
end
if ~isfield(opt, 'max_it') || isempty(opt.max_it)
    opt.max_it = 150;
end
if ~isfield(opt, 'sc') || ~isfield(opt.sc, 'red_it') || isempty(opt.sc.red_it)
    opt.sc.red_it = 20;
end
if ~isfield(opt, 'step_control') || isempty(opt.step_control)
    opt.step_control = 0;
end
if ~isfield(opt, 'cost_mult') || isempty(opt.cost_mult)
    opt.cost_mult = 1;
end
if ~isfield(opt, 'verbose') || isempty(opt.verbose)
    opt.verbose = 0;
end
%% algorithm constants
if ~isfield(opt, 'xi') || isempty(opt.xi)
    opt.xi = 0.99995;           %% OPT_IPM_PHI
end
if ~isfield(opt, 'sigma') || isempty(opt.sigma)
    opt.sigma = 0.1;            %% OPT_IPM_SIGMA, centering parameter
end
if ~isfield(opt, 'z0') || isempty(opt.z0)
    opt.z0 = 1;                 %% OPT_IPM_INIT_SLACK
end
if ~isfield(opt, 'alpha_min') || isempty(opt.alpha_min)
    opt.alpha_min = 1e-8;       %% OPT_AP_AD_MIN
end
if ~isfield(opt, 'rho_min') || isempty(opt.rho_min)
    opt.rho_min = 0.95;         %% OPT_IPM_QUAD_LOWTHRESH
end
if ~isfield(opt, 'rho_max') || isempty(opt.rho_max)
    opt.rho_max = 1.05;         %% OPT_IPM_QUAD_HIGHTHRESH
end
if ~isfield(opt, 'mu_threshold') || isempty(opt.mu_threshold)
    opt.mu_threshold = 1e-5;    %% SCOPF_MULTIPLIERS_FILTER_THRESH
end
if ~isfield(opt, 'max_stepsize') || isempty(opt.max_stepsize)
    opt.max_stepsize = 1e10;
end

%% initialize history
hist(opt.max_it+1) = struct('feascond', 0, 'gradcond', 0, 'compcond', 0, ...
    'costcond', 0, 'gamma', 0, 'stepsize', 0, 'obj', 0, ...
    'alphap', 0, 'alphad', 0);

%%-----  set up problem  -----
%% constants
[xi, sigma, z0, alpha_min, rho_min, rho_max, mu_threshold, max_stepsize] = ...
    deal(opt.xi, opt.sigma, opt.z0, opt.alpha_min, ...
        opt.rho_min, opt.rho_max, opt.mu_threshold, opt.max_stepsize);
if xi >= 1 || xi < 0.5
    error('mips: opt.xi (%g) must be a number slightly less than 1', opt.xi);
end
if sigma > 1 || sigma <= 0
    error('mips: opt.sigma (%g) must be a number between 0 and 1', opt.sigma);
end

%% initialize
i = 0;                      %% iteration counter
converged = 0;              %% flag
eflag = 0;                  %% exit flag

%% add var limits to linear constraints
AA = [speye(nx); A];
ll = [xmin; l];
uu = [xmax; u];

%% split up linear constraints
ieq = find( abs(uu-ll) <= eps );            %% equality
igt = find( uu >=  1e10 & ll > -1e10 );     %% greater than, unbounded above
ilt = find( ll <= -1e10 & uu <  1e10 );     %% less than, unbounded below
ibx = find( (abs(uu-ll) > eps) & (uu < 1e10) & (ll > -1e10) );
Ae = AA(ieq, :);
be = uu(ieq, 1);
Ai  = [ AA(ilt, :); -AA(igt, :); AA(ibx, :); -AA(ibx, :) ];
bi  = [ uu(ilt, 1); -ll(igt, 1); uu(ibx, 1); -ll(ibx, 1) ];

%% evaluate cost f(x0) and constraints g(x0), h(x0)
x = x0;
[f, df] = f_fcn(x);             %% cost
f = f * opt.cost_mult;
df = df * opt.cost_mult;
if nonlinear
    [hn, gn, dhn, dgn] = gh_fcn(x); %% nonlinear constraints
    h = [hn; Ai * x - bi];          %% inequality constraints
    g = [gn; Ae * x - be];          %% equality constraints
    dh = [dhn Ai'];                 %% 1st derivative of inequalities
    dg = [dgn Ae'];                 %% 1st derivative of equalities
else
    h = Ai * x - bi;                %% inequality constraints
    g = Ae * x - be;                %% equality constraints
    dh = Ai';                       %% 1st derivative of inequalities
    dg = Ae';                       %% 1st derivative of equalities
end

%% grab some dimensions
neq = size(g, 1);           %% number of equality constraints
niq = size(h, 1);           %% number of inequality constraints
neqnln = size(gn, 1);       %% number of nonlinear equality constraints
niqnln = size(hn, 1);       %% number of nonlinear inequality constraints
nlt = length(ilt);          %% number of upper bounded linear inequalities
ngt = length(igt);          %% number of lower bounded linear inequalities
nbx = length(ibx);          %% number of doubly bounded linear inequalities

%% initialize gamma, lam, mu, z, e
gamma = 1;                  %% barrier coefficient, r in Harry's code
lam = zeros(neq, 1);
z   = z0 * ones(niq, 1);
mu  = z;
k = find(h < -z0);
z(k) = -h(k);
k = find(gamma / z > z0);   %% (seems k is always empty if gamma = z0 = 1)
if ~isempty(k)
    mu(k) = gamma / z(k);
end
e = ones(niq, 1);

%% check tolerance
f0 = f;
if opt.step_control
    L = f + lam' * g + mu' * (h+z) - gamma * sum(log(z));
end
Lx = df + dg * lam + dh * mu;
if isempty(h)
    maxh = zeros(1,0);
else
    maxh = max(h);
end
feascond = max([norm(g, Inf), maxh]) / (1 + max([ norm(x, Inf), norm(z, Inf) ]));
gradcond = norm(Lx, Inf) / (1 + max([ norm(lam, Inf), norm(mu, Inf) ]));
compcond = (z' * mu) / (1 + norm(x, Inf));
costcond = abs(f - f0) / (1 + abs(f0));
%% save history
hist(i+1) = struct('feascond', feascond, 'gradcond', gradcond, ...
    'compcond', compcond, 'costcond', costcond, 'gamma', gamma, ...
    'stepsize', 0, 'obj', f/opt.cost_mult, 'alphap', 0, 'alphad', 0);
if strcmp(upper(opt.linsolver), 'PARDISO')
    ls = 'PARDISO';
    mplinsolve_opt = struct('pardiso', ...
                            struct('mtype', -2, ...
                                   'iparm', [8, 1;  %% max it refinement steps
                                             10, 12;%% eps pivot
                                             11, 1; %% use scaling vectors
                                             13, 2; %% improved accuracy
                                             18, 0; %% do not determine nnz in LU
                                             21, 3; %% ? undocumented pivoting
                                             ] ));
    if ~have_fcn('pardiso')
        warning('mips: PARDISO linear solver not available, using default');
        opt.linsolver = '';
        ls = 'built-in';
    end
else
    ls = 'built-in';
    mplinsolve_opt = [];
end
if opt.verbose
    if opt.step_control, s = '-sc'; else, s = ''; end
    v = mipsver('all');
    fprintf('MATPOWER Interior Point Solver -- MIPS%s, Version %s, %s\n (using %s linear solver)', ...
        s, v.Version, v.Date, ls);
    if opt.verbose > 1
        fprintf('\n it    objective   step size   feascond     gradcond     compcond     costcond  ');
        fprintf('\n----  ------------ --------- ------------ ------------ ------------ ------------');
        fprintf('\n%3d  %12.8g %10s %12g %12g %12g %12g', ...
            i, f/opt.cost_mult, '', feascond, gradcond, compcond, costcond);
    end
end
if feascond < opt.feastol && gradcond < opt.gradtol && ...
                compcond < opt.comptol && costcond < opt.costtol
    converged = 1;
    if opt.verbose
        fprintf('\nConverged!\n');
    end
end

%%-----  do Newton iterations  -----
while (~converged && i < opt.max_it)
    %% update iteration counter
    i = i + 1;

    %% compute update step
    lambda = struct('eqnonlin', lam(1:neqnln), 'ineqnonlin', mu(1:niqnln));
    if nonlinear
        if isempty(hess_fcn)
            fprintf('mips: Hessian evaluation via finite differences not yet implemented.\n       Please provide your own hessian evaluation function.');
        end
        Lxx = hess_fcn(x, lambda, opt.cost_mult);
    else
        [f_, df_, d2f] = f_fcn(x);      %% cost
        Lxx = d2f * opt.cost_mult;
    end
    zinvdiag = sparse(1:niq, 1:niq, 1 ./ z, niq, niq);
    mudiag = sparse(1:niq, 1:niq, mu, niq, niq);
    dh_zinv = dh * zinvdiag;
    M = Lxx + dh_zinv * mudiag * dh';
    N = Lx + dh_zinv * (mudiag * h + gamma * e);
    dxdlam = mplinsolve([M dg; dg' sparse(neq, neq)], [-N; -g], opt.linsolver, mplinsolve_opt);
%     AAA = [
%         M  dg;
%         dg'  sparse(neq, neq)
%     ];
%     rc = 1/condest(AAA);
%     if rc < 1e-22
%         fprintf('my RCOND = %g\n', rc);
%         n = size(AAA, 1);
%         AAA = AAA + 1e-3 * speye(n,n);
%     end
%     bbb = [-N; -g];
%     dxdlam = AAA \ bbb;
    if any(isnan(dxdlam)) || norm(dxdlam) > max_stepsize
        if opt.verbose
            fprintf('\nNumerically Failed\n');
        end
        eflag = -1;
        break;
    end
    dx = dxdlam(1:nx, 1);
    dlam = dxdlam(nx+(1:neq), 1);
    dz = -h - z - dh' * dx;
    dmu = -mu + zinvdiag *(gamma*e - mudiag * dz);

    %% optional step-size control
    sc = 0;
    if opt.step_control
        x1 = x + dx;

        %% evaluate cost, constraints, derivatives at x1
        [f1, df1] = f_fcn(x1);          %% cost
        f1 = f1 * opt.cost_mult;
        df1 = df1 * opt.cost_mult;
        if nonlinear
            [hn1, gn1, dhn1, dgn1] = gh_fcn(x1); %% nonlinear constraints
            h1 = [hn1; Ai * x1 - bi];       %% inequality constraints
            g1 = [gn1; Ae * x1 - be];       %% equality constraints
            dh1 = [dhn1 Ai'];               %% 1st derivative of inequalities
            dg1 = [dgn1 Ae'];               %% 1st derivative of equalities
        else
            h1 = Ai * x1 - bi;              %% inequality constraints
            g1 = Ae * x1 - be;              %% equality constraints
            dh1 = dh;                       %% 1st derivative of inequalities
            dg1 = dg;                       %% 1st derivative of equalities
        end

        %% check tolerance
        Lx1 = df1 + dg1 * lam + dh1 * mu;
        if isempty(h1)
            maxh1 = zeros(1,0);
        else
            maxh1 = max(h1);
        end
        feascond1 = max([norm(g1, Inf), maxh1]) / (1 + max([ norm(x1, Inf), norm(z, Inf) ]));
        gradcond1 = norm(Lx1, Inf) / (1 + max([ norm(lam, Inf), norm(mu, Inf) ]));

        if feascond1 > feascond && gradcond1 > gradcond
            sc = 1;
        end
    end
    if sc
        alpha = 1;
        for j = 1:opt.sc.red_it
            dx1 = alpha * dx;
            x1 = x + dx1;
            f1 = f_fcn(x1);                 %% cost
            f1 = f1 * opt.cost_mult;
            if nonlinear
                [hn1, gn1] = gh_fcn(x1);    %% nonlinear constraints
                h1 = [hn1; Ai * x1 - bi];   %% inequality constraints
                g1 = [gn1; Ae * x1 - be];   %% equality constraints
            else
                h1 = Ai * x1 - bi;          %% inequality constraints
                g1 = Ae * x1 - be;          %% equality constraints
            end
            L1 = f1 + lam' * g1 + mu' * (h1+z) - gamma * sum(log(z));
            if opt.verbose > 2
                fprintf('\n   %3d            %10g', -j, norm(dx1));
            end
            rho = (L1 - L) / (Lx' * dx1 + 0.5 * dx1' * Lxx * dx1);
            if rho > rho_min && rho < rho_max
                break;
            else
                alpha = alpha / 2;
            end
        end
        dx = alpha * dx;
        dz = alpha * dz;
        dlam = alpha * dlam;
        dmu = alpha * dmu;
    end

    %% do the update
    k = find(dz < 0);
    if isempty(k)
        alphap = 1;
    else
        alphap = min( [xi * min(z(k) ./ -dz(k)) 1] );
    end
    k = find(dmu < 0);
    if isempty(k)
        alphad = 1;
    else
        alphad = min( [xi * min(mu(k) ./ -dmu(k)) 1] );
    end
    x = x + alphap * dx;
    z = z + alphap * dz;
    lam = lam + alphad * dlam;
    mu  = mu  + alphad * dmu;
    if niq > 0
        gamma = sigma * (z' * mu) / niq;
    end

    %% evaluate cost, constraints, derivatives
    [f, df] = f_fcn(x);                 %% cost
    f = f * opt.cost_mult;
    df = df * opt.cost_mult;
    if nonlinear
        [hn, gn, dhn, dgn] = gh_fcn(x); %% nonlinear constraints
        h = [hn; Ai * x - bi];          %% inequality constraints
        g = [gn; Ae * x - be];          %% equality constraints
        dh = [dhn Ai'];                 %% 1st derivative of inequalities
        dg = [dgn Ae'];                 %% 1st derivative of equalities
    else
        h = Ai * x - bi;                %% inequality constraints
        g = Ae * x - be;                %% equality constraints
        %% 1st derivatives are constant, still dh = Ai', dg = Ae'
    end

    %% check tolerance
    Lx = df + dg * lam + dh * mu;
    if isempty(h)
        maxh = zeros(1,0);
    else
        maxh = max(h);
    end
    feascond = max([norm(g, Inf), maxh]) / (1 + max([ norm(x, Inf), norm(z, Inf) ]));
    gradcond = norm(Lx, Inf) / (1 + max([ norm(lam, Inf), norm(mu, Inf) ]));
    compcond = (z' * mu) / (1 + norm(x, Inf));
    costcond = abs(f - f0) / (1 + abs(f0));
    %% save history
    hist(i+1) = struct('feascond', feascond, 'gradcond', gradcond, ...
        'compcond', compcond, 'costcond', costcond, 'gamma', gamma, ...
        'stepsize', norm(dx), 'obj', f/opt.cost_mult, ...
        'alphap', alphap, 'alphad', alphad);

    if opt.verbose > 1
        fprintf('\n%3d  %12.8g %10.5g %12g %12g %12g %12g', ...
            i, f/opt.cost_mult, norm(dx), feascond, gradcond, compcond, costcond);
    end
    if feascond < opt.feastol && gradcond < opt.gradtol && ...
                    compcond < opt.comptol && costcond < opt.costtol
        converged = 1;
        if opt.verbose
            fprintf('\nConverged!\n');
        end
    else
        if any(isnan(x)) || alphap < alpha_min || alphad < alpha_min || ...
                gamma < eps || gamma > 1/eps
            if opt.verbose
                fprintf('\nNumerically Failed\n');
            end
            eflag = -1;
            break;
        end
        f0 = f;
        if opt.step_control
            L = f + lam' * g + mu' * (h+z) - gamma * sum(log(z));
        end
    end
end

if opt.verbose
    if ~converged
        fprintf('\nDid not converge in %d iterations.\n', i);
    end
end

%%-----  package up results  -----
hist = hist(1:i+1);
if eflag ~= -1
    eflag = converged;
end
output = struct('iterations', i, 'hist', hist, 'message', '');
if eflag == 0
    output.message = 'Did not converge';
elseif eflag == 1
    output.message = 'Converged';
elseif eflag == -1
    output.message = 'Numerically failed';
else
    output.message = 'Please hang up and dial again';
end

%% zero out multipliers on non-binding constraints
mu(h < -opt.feastol & mu < mu_threshold) = 0;

%% un-scale cost and prices
f   = f   / opt.cost_mult;
lam = lam / opt.cost_mult;
mu  = mu  / opt.cost_mult;

%% re-package multipliers into struct
lam_lin = lam((neqnln+1):neq);              %% lambda for linear constraints
mu_lin  = mu((niqnln+1):niq);               %% mu for linear constraints
kl = find(lam_lin < 0);                     %% lower bound binding
ku = find(lam_lin > 0);                     %% upper bound binding

mu_l = zeros(nx+nA, 1);
mu_l(ieq(kl)) = -lam_lin(kl);
mu_l(igt) = mu_lin(nlt+(1:ngt));
mu_l(ibx) = mu_lin(nlt+ngt+nbx+(1:nbx));

mu_u = zeros(nx+nA, 1);
mu_u(ieq(ku)) = lam_lin(ku);
mu_u(ilt) = mu_lin(1:nlt);
mu_u(ibx) = mu_lin(nlt+ngt+(1:nbx));

fields = { ...
    'mu_l', mu_l((nx+1):end), ...
    'mu_u', mu_u((nx+1):end), ...
    'lower', mu_l(1:nx), ...
    'upper', mu_u(1:nx) ...
};

if niqnln > 0
    fields = { ...
        'ineqnonlin', mu(1:niqnln), ...
        fields{:} ...
    };
end
if neqnln > 0
    fields = { ...
        'eqnonlin', lam(1:neqnln), ...
        fields{:} ...
    };
end

lambda = struct(fields{:});

% lambda = struct( ...
%     'eqnonlin', lam(1:neqnln), ...
%     'ineqnonlin', mu(1:niqnln), ...
%     'mu_l', mu_l((nx+1):end), ...
%     'mu_u', mu_u((nx+1):end), ...
%     'lower', mu_l(1:nx), ...
%     'upper', mu_u(1:nx) );
