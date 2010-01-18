function [x, f, eflag, output, lambda] = pdipm(ipm_f, x0, A, l, u, xmin, xmax, ipm_gh, ipm_hess, opt)
%PDIPM  Primal-dual interior point method for NLP.
%   Minimize a function f beginning from a starting point x0, subject to
%   optional linear and non-linear constraints and variable bounds.
%
%       min f(x)
%        x
%
%   such that
%
%           h(x) = 0            (non-linear equalities)
%           g(x) <= 0           (non-linear inequalities)
%           l <= A*x <= u       (linear constraints)
%           xmin <= x <= xmax   (variable bounds)
%
%   [x, f, exitflag, output, lambda] = ...
%       pdipm(f, x0, A, l, u, xmin, xmax, gh, hess, opt)
%   x = pdipm(f, x0)
%   x = pdipm(f, x0, A, l)
%   x = pdipm(f, x0, A, l, u)
%   x = pdipm(f, x0, A, l, u, xmin)
%   x = pdipm(f, x0, A, l, u, xmin, xmax)
%   x = pdipm(f, x0, A, l, u, xmin, xmax, gh)
%   x = pdipm(f, x0, A, l, u, xmin, xmax, gh, hess)
%   x = pdipm(f, x0, A, l, u, xmin, xmax, gh, hess, opt)
%   x = pdipm(problem), where problem is a struct with fields:
%                       f, x0, A, l, u, xmin, xmax, gh, hess, opt
%                       all fields except 'f' and 'x0' are optional
%   x = pdipm(...)
%   [x, fval] = pdipm(...)
%   [x, fval, exitflag] = pdipm(...)
%   [x, fval, exitflag, output] = pdipm(...)
%   [x, fval, exitflag, output, lambda] = pdipm(...)
%
%   Inputs:
%       f : handle to function that evaluates the objective function,
%           it's gradients and Hessian for a given value of x. If there
%           are non-linear constraints the Hessian information is
%           provided by the 'hess' function passed in the 9th argument
%           and is not required here. Calling syntax for this function:
%               [f, df, d2f] = f(x)
%       x0 : starting value of optimization vector x
%       A, l, u : define the optional linear constraints. Default
%           values for the elements of l and u are -Inf and Inf,
%           respectively.
%       xmin, xmax : optional lower and upper bounds on the
%           x variables, defaults are -Inf and Inf, respectively.
%       gh : handle to function that evaluates the optional
%           non-linear constraints and their gradients for a given
%           value of x. Calling syntax for this function is:
%               [g, h, dg, dh] = gh(x)
%       hess : handle to function that computes the Hessian of the
%           Lagrangian for given values of x, lambda and mu, where
%           lambda and mu are the multipliers on the equality and
%           inequality constraints, h and g, respectively. The calling
%           syntax for this function is:
%               Lxx = hess(x, lam)
%           where lambda = lam.eqnonlin and mu = lam.ineqnonlin.
%       opt : optional options structure with the following fields,
%           all of which are also optional (default values in shown in
%           parentheses)
%           verbose (0) - controls level of progress output displayed
%           feastol (1e-6) - termination tolerance for feasibility
%               condition
%           gradtol (1e-6) - termination tolerance for gradient
%               condition
%           comptol (1e-6) - termination tolerance for complementarity
%               condition
%           costtol (1e-6) - termination tolerance for cost condition
%           max_it (150) - maximum number of iterations
%           step_control (0) - set to 1 to turn on step-size control
%           max_red (20) - maximum number of step-size reductions if
%               step-control is on
%           cost_mult (1) - cost multiplier used to scale the objective
%               function for improved conditioning. Note: The same
%               value must also be passed to the Hessian evaluation
%               function so that it can appropriately scale the
%               objective function term in the Hessian of the
%               Lagrangian.
%       problem : The inputs can alternatively be supplied in a single
%           struct with fields corresponding to the input arguments
%           described above: f, x0, A, l, u, xmin, xmax, gh, hess, opt
%
%   Outputs:
%       x : solution vector
%       fval : final objective function value
%       exitflag : exit flag,
%           1 = first order optimality conditions satisfied
%           0 = maximum number of iterations reached
%           -1 = numerically failed
%       output : structure with fields:
%           iterations - number of iterations performed
%           hist - struct array with trajectories of the following:
%                   feascond, gradcond, compcond, costcond, gamma,
%                   stepsize, obj, alphap, alphad
%       lambda : struct containing the Langrange and Kuhn-Tucker
%           multipliers on the constraints, with fields:
%           eqnonlin - non-linear equality constraints
%           ineqnonlin - non-linear inequality constraints
%           mu_l - lower bound on linear constraints
%           mu_u - upper bound on linear constraints
%           lower - lower bound on optimization variables
%           upper - upper bound on optimization variables
%
%   Note the calling syntax is almost identical to that of 'fmincon'
%   from MathWorks' Optimization Toolbox. The main difference is that
%   the linear constraints are specified with A, l, u instead of
%   A, b, Aeq, beq. The functions for evaluating the objective
%   function, constraints and Hessian are identical.
%
%   Ported by Ray Zimmerman from C code written by H. Wang for his
%   PhD dissertation:
%     "On the Computation and Application of Multi-period
%     Security-Constrained Optimal Power Flow for Real-time
%     Electricity Market Operations", Cornell University, May 2007.
%
%   See also:
%     H. Wang, C. E. Murillo-S‡nchez, R. D. Zimmerman, R. J. Thomas,
%     "On Computational Issues of Market-Based Optimal Power Flow",
%     IEEE Transactions on Power Systems, Vol. 22, No. 3, Aug. 2007,
%     pp. 1185-1193.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009-2010 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- input argument handling  -----
%% gather inputs
if nargin == 1 && isstruct(ipm_f)       %% problem struct
    p = ipm_f;
    ipm_f = p.f;
    x0 = p.x0;
    nx = size(x0, 1);       %% number of optimization variables
    if isfield(p, 'opt'),   opt = p.opt;        else,   opt = [];       end
    if isfield(p, 'hess'),  ipm_hess = p.hess;  else,   ipm_hess = '';  end
    if isfield(p, 'gh'),    ipm_gh = p.gh;      else,   ipm_gh = '';    end
    if isfield(p, 'xmax'),  xmax = p.xmax;      else,   xmax = [];      end
    if isfield(p, 'xmin'),  xmin = p.xmin;      else,   xmin = [];      end
    if isfield(p, 'u'),     u = p.u;            else,   u = [];         end
    if isfield(p, 'l'),     l = p.l;            else,   l = [];         end
    if isfield(p, 'A'),     A = p.A;            else,   A=sparse(0,nx); end
else                                    %% individual args
    nx = size(x0, 1);       %% number of optimization variables
    if nargin < 10
        opt = [];
        if nargin < 9
            ipm_hess = '';
            if nargin < 8
                ipm_gh = '';
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
if isempty(ipm_gh)
    nonlinear = false;          %% no non-linear constraints present
    gn = []; hn = [];
else
    nonlinear = true;           %% we have some non-linear constraints
end

%% default options
if isempty(opt)
    opt = struct('verbose', 0);
end
if ~isfield(opt, 'feastol') || isempty(opt.feastol)
    opt.feastol = 1e-6;
end
if ~isfield(opt, 'gradtol') || isempty(opt.gradtol)
    opt.gradtol = 1e-6;
end
if ~isfield(opt, 'comptol') || isempty(opt.comptol)
    opt.comptol = 1e-6;
end
if ~isfield(opt, 'costtol') || isempty(opt.costtol)
    opt.costtol = 1e-6;
end
if ~isfield(opt, 'max_it') || isempty(opt.max_it)
    opt.max_it = 150;
end
if ~isfield(opt, 'max_red') || isempty(opt.max_red)
    opt.max_red = 20;
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

%% initialize history
hist(opt.max_it+1) = struct('feascond', 0, 'gradcond', 0, 'compcond', 0, ...
    'costcond', 0, 'gamma', 0, 'stepsize', 0, 'obj', 0, ...
    'alphap', 0, 'alphad', 0);

%%-----  set up problem  -----
%% constants
xi = 0.99995;           %% OPT_IPM_PHI
sigma = 0.1;            %% OPT_IPM_SIGMA
z0 = 1;                 %% OPT_IPM_INIT_SLACK
alpha_min = 1e-8;       %% OPT_AP_AD_MIN
rho_min = 0.95;         %% OPT_IPM_QUAD_LOWTHRESH
rho_max = 1.05;         %% OPT_IPM_QUAD_HIGHTHRESH
mu_threshold = 1e-5;    %% SCOPF_MULTIPLIERS_FILTER_THRESH

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
be = uu(ieq);
Ai  = [ AA(ilt, :); -AA(igt, :); AA(ibx, :); -AA(ibx, :) ];
bi  = [ uu(ilt);   -ll(igt);     uu(ibx);    -ll(ibx)];

%% evaluate cost f(x0) and constraints g(x0), h(x0)
x = x0;
[f, df] = ipm_f(x);             %% cost
f = f * opt.cost_mult;
df = df * opt.cost_mult;
if nonlinear
    [gn, hn, dgn, dhn] = ipm_gh(x); %% non-linear constraints
    g = [gn; Ai * x - bi];          %% inequality constraints
    h = [hn; Ae * x - be];          %% equality constraints
    dg = [dgn Ai'];                 %% 1st derivative of inequalities
    dh = [dhn Ae'];                 %% 1st derivative of equalities
else
    g = Ai * x - bi;                %% inequality constraints
    h = Ae * x - be;                %% equality constraints
    dg = Ai';                       %% 1st derivative of inequalities
    dh = Ae';                       %% 1st derivative of equalities
end

%% grab some dimensions
neq = size(h, 1);           %% number of equality constraints
niq = size(g, 1);           %% number of inequality constraints
neqnln = size(hn, 1);       %% number of non-linear equality constraints
niqnln = size(gn, 1);       %% number of non-linear inequality constraints
nlt = length(ilt);          %% number of upper bounded linear inequalities
ngt = length(igt);          %% number of lower bounded linear inequalities
nbx = length(ibx);          %% number of doubly bounded linear inequalities

%% initialize gamma, lam, mu, z, e
gamma = 1;                  %% barrier coefficient, r in Harry's code
lam = zeros(neq, 1);
z   = z0 * ones(niq, 1);
mu  = z;
k = find(g < -z0);
z(k) = -g(k);
k = find(gamma / z > z0);   %% (seems k is always empty if gamma = z0 = 1)
if ~isempty(k)
    mu(k) = gamma / z(k);
end
e = ones(niq, 1);

%% check tolerance
f0 = f;
if opt.step_control
    L = f + lam' * h + mu' * (g+z) - gamma * sum(log(z));
end
Lx = df + dh * lam + dg * mu;
feascond = max([norm(h, Inf), max(g)]) / (1 + max([ norm(x, Inf), norm(z, Inf) ]));
gradcond = norm(Lx, Inf) / (1 + max([ norm(lam, Inf), norm(mu, Inf) ]));
compcond = (z' * mu) / (1 + norm(x, Inf));
costcond = abs(f - f0) / (1 + abs(f0));
%% save history
hist(i+1) = struct('feascond', feascond, 'gradcond', gradcond, ...
    'compcond', compcond, 'costcond', costcond, 'gamma', gamma, ...
    'stepsize', 0, 'obj', f/opt.cost_mult, 'alphap', 0, 'alphad', 0);
if opt.verbose > 1
    fprintf('\n it    objective   step size   feascond     gradcond     compcond     costcond  ');
    fprintf('\n----  ------------ --------- ------------ ------------ ------------ ------------');
    fprintf('\n%3d  %12.8g %10s %12g %12g %12g %12g', ...
        i, f/opt.cost_mult, '', feascond, gradcond, compcond, costcond);
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
        if isempty(ipm_hess)
            fprintf('pdipm: Hessian evaluation via finite differences not yet implemented.\n       Please provide your own hessian evaluation function.');
        end
        Lxx = ipm_hess(x, lambda);
    else
        [f_, df_, d2f] = ipm_f(x);      %% cost
        Lxx = d2f * opt.cost_mult;
    end
    zinvdiag = sparse(1:niq, 1:niq, 1 ./ z, niq, niq);
    mudiag = sparse(1:niq, 1:niq, mu, niq, niq);
    dg_zinv = dg * zinvdiag;
    M = Lxx + dg_zinv * mudiag * dg';
    N = Lx + dg_zinv * (mudiag * g + gamma * e);
    dxdlam = [M dh; dh' sparse(neq, neq)] \ [-N; -h];
%     AAA = [
%         M  dh;
%         dh'  sparse(neq, neq)
%     ];
%     rc = 1/condest(AAA);
%     if rc < 1e-22
%         fprintf('my RCOND = %g\n', rc);
%         n = size(AAA, 1);
%         AAA = AAA + 1e-3 * speye(n,n);
%     end
%     bbb = [-N; -h];
%     dxdlam = AAA \ bbb;
    dx = dxdlam(1:nx);
    dlam = dxdlam(nx+(1:neq));
    dz = -g - z - dg' * dx;
    dmu = -mu + zinvdiag *(gamma*e - mudiag * dz);

    %% optional step-size control
    sc = 0;
    if opt.step_control
        x1 = x + dx;

        %% evaluate cost, constraints, derivatives at x1
        [f1, df1] = ipm_f(x1);          %% cost
        f1 = f1 * opt.cost_mult;
        df1 = df1 * opt.cost_mult;
        if nonlinear
            [gn1, hn1, dgn1, dhn1] = ipm_gh(x1); %% non-linear constraints
            g1 = [gn1; Ai * x1 - bi];       %% inequality constraints
            h1 = [hn1; Ae * x1 - be];       %% equality constraints
            dg1 = [dgn1 Ai'];               %% 1st derivative of inequalities
            dh1 = [dhn1 Ae'];               %% 1st derivative of equalities
        else
            g1 = Ai * x1 - bi;              %% inequality constraints
            h1 = Ae * x1 - be;              %% equality constraints
            dg1 = dg;                       %% 1st derivative of inequalities
            dh1 = dh;                       %% 1st derivative of equalities
        end

        %% check tolerance
        Lx1 = df1 + dh1 * lam + dg1 * mu;
        feascond1 = max([norm(h1, Inf), max(g1)]) / (1 + max([ norm(x1, Inf), norm(z, Inf) ]));
        gradcond1 = norm(Lx1, Inf) / (1 + max([ norm(lam, Inf), norm(mu, Inf) ]));

        if feascond1 > feascond && gradcond1 > gradcond
            sc = 1;
        end
    end
    if sc
        alpha = 1;
        for j = 1:opt.max_red
            dx1 = alpha * dx;
            x1 = x + dx1;
            f1 = ipm_f(x1);             %% cost
            f1 = f1 * opt.cost_mult;
            if nonlinear
                [gn1, hn1] = ipm_gh(x1);    %% non-linear constraints
                g1 = [gn1; Ai * x1 - bi];   %% inequality constraints
                h1 = [hn1; Ae * x1 - be];   %% equality constraints
            else
                g1 = Ai * x1 - bi;          %% inequality constraints
                h1 = Ae * x1 - be;          %% equality constraints
            end
            L1 = f1 + lam' * h1 + mu' * (g1+z) - gamma * sum(log(z));
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
    alphap = min( [xi * min(z(k) ./ -dz(k)) 1] );
    k = find(dmu < 0);
    alphad = min( [xi * min(mu(k) ./ -dmu(k)) 1] );
    x = x + alphap * dx;
    z = z + alphap * dz;
    lam = lam + alphad * dlam;
    mu  = mu  + alphad * dmu;
    gamma = sigma * (z' * mu) / niq;

    %% evaluate cost, constraints, derivatives
    [f, df] = ipm_f(x);             %% cost
    f = f * opt.cost_mult;
    df = df * opt.cost_mult;
    if nonlinear
        [gn, hn, dgn, dhn] = ipm_gh(x); %% non-linear constraints
        g = [gn; Ai * x - bi];          %% inequality constraints
        h = [hn; Ae * x - be];          %% equality constraints
        dg = [dgn Ai'];                 %% 1st derivative of inequalities
        dh = [dhn Ae'];                 %% 1st derivative of equalities
    else
        g = Ai * x - bi;                %% inequality constraints
        h = Ae * x - be;                %% equality constraints
        %% 1st derivatives are constant, still dg = Ai', dh = Ae'
    end

    %% check tolerance
    Lx = df + dh * lam + dg * mu;
    feascond = max([norm(h, Inf), max(g)]) / (1 + max([ norm(x, Inf), norm(z, Inf) ]));
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
            L = f + lam' * h + mu' * (g+z) - gamma * sum(log(z));
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
output = struct('iterations', i, 'hist', hist);

%% zero out multipliers on non-binding constraints
mu(g < -opt.feastol & mu < mu_threshold) = 0;

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

fields = {};
if neqnln > 0
    fields = {fields{:}, 'eqnonlin', lam(1:neqnln)};
end
if niqnln > 0
    fields = {fields{:}, 'ineqnonlin', mu(1:niqnln)};
end
if neq > neqnln || niq > niqnln
    fields = {fields{:}, 'mu_l', mu_l((nx+1):end), 'mu_u', mu_u((nx+1):end)};
end
if any(xmin ~= -Inf)
    fields = {fields{:}, 'lower', mu_l(1:nx)};
end
if any(xmax ~= Inf)
    fields = {fields{:}, 'upper', mu_u(1:nx)};
end

lambda = struct(fields{:});

% lambda = struct( ...
%     'eqnonlin', lam(1:neqnln), ...
%     'ineqnonlin', mu(1:niqnln), ...
%     'mu_l', mu_l((nx+1):end), ...
%     'mu_u', mu_u((nx+1):end), ...
%     'lower', mu_l(1:nx), ...
%     'upper', mu_u(1:nx) );
