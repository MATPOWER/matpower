function [x, f, eflag, output, lambda] = nlps_ipopt(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt)
% nlps_ipopt - Nonlinear programming (NLP) Solver based on IPOPT.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       NLPS_IPOPT(F_FCN, X0, A, L, U, XMIN, XMAX, GH_FCN, HESS_FCN, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = NLPS_IPOPT(PROBLEM)
%   A wrapper function providing a standardized interface for using
%   IPOPT to solve the following NLP (nonlinear programming) problem:
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
%           ipopt_opt   - options struct for IPOPT, value in verbose
%                   overrides these options
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
%           0 = failed to converge
%       OUTPUT : output struct with the following fields:
%           status - see IPOPT documentation for INFO.status
%             https://coin-or.github.io/Ipopt/IpReturnCodes__inc_8h_source.html
%           iterations - number of iterations performed (INFO.iter)
%           cpu - see IPOPT documentation for INFO.cpu
%           eval - see IPOPT documentation for INFO.eval
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
%           nlps_ipopt(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
%
%       x = nlps_ipopt(f_fcn, x0);
%       x = nlps_ipopt(f_fcn, x0, A, l);
%       x = nlps_ipopt(f_fcn, x0, A, l, u);
%       x = nlps_ipopt(f_fcn, x0, A, l, u, xmin);
%       x = nlps_ipopt(f_fcn, x0, A, l, u, xmin, xmax);
%       x = nlps_ipopt(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn);
%       x = nlps_ipopt(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn);
%       x = nlps_ipopt(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
%       x = nlps_ipopt(problem);
%               where problem is a struct with fields:
%                   f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt
%                   all fields except 'f_fcn' and 'x0' are optional
%       x = nlps_ipopt(...);
%       [x, f] = nlps_ipopt(...);
%       [x, f, exitflag] = nlps_ipopt(...);
%       [x, f, exitflag, output] = nlps_ipopt(...);
%       [x, f, exitflag, output, lambda] = nlps_ipopt(...);
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
%       [x, f, exitflag, output, lambda] = nlps_ipopt(problem);
%
% See also nlps_master, ipopt.

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

%% make sure A is sparse
if ~issparse(A)
    A = sparse(A);
end

%% replace equality variable bounds with an equality constraint
%% (since IPOPT does not return shadow prices on variables that it eliminates)
kk = find(xmin == xmax);
nk = length(kk);
if nk
    A = [ A; sparse((1:nk)', kk, 1, nk, nx) ];
    l = [ l; xmin(kk) ];
    u = [ u; xmax(kk) ];
    xmin(kk) = -Inf;
    xmax(kk) = Inf;
    nA = size(A, 1);            %% updated number of linear constraints
end

%%-----  set up problem  -----
%% build Jacobian and Hessian structure
randx = rand(size(x0));
nonz = 1e-20;
if nonlinear
    [h, g, dhs, dgs] = gh_fcn(randx);
    dhs(dhs ~= 0) = nonz;   %% set non-zero entries to tiny value (for adding later)
    dgs(dgs ~= 0) = nonz;   %% set non-zero entries to tiny value (for adding later)
    Js = [dgs'; dhs'; A];
else
    g = []; h = [];
    dhs = sparse(0, nx); dgs = dhs; Js = A;
end
neq = length(g);
niq = length(h);
if isempty(hess_fcn)
    [f_, df_, Hs] = f_fcn(randx);   %% cost
else
    lam = struct('eqnonlin', rand(neq, 1), 'ineqnonlin', rand(niq, 1));
    Hs = hess_fcn(randx, lam, 1);
end
if ~issparse(Hs), Hs = sparse(Hs); end      %% convert to sparse if necessary
Hs(Hs ~= 0) = nonz;     %% set non-zero entries to tiny value (for adding later)

%% set options struct for IPOPT
if ~isempty(opt) && isfield(opt, 'ipopt_opt') && ~isempty(opt.ipopt_opt)
    options.ipopt = ipopt_options(opt.ipopt_opt);
else
    options.ipopt = ipopt_options;
end
if verbose
    options.ipopt.print_level = min(12, verbose*2+1);
else
    options.ipopt.print_level = 0;
end

%% extra data to pass to functions
options.auxdata = struct( ...
    'f_fcn',    f_fcn, ...
    'gh_fcn',   gh_fcn, ...
    'hess_fcn', hess_fcn, ...
    'A',        A, ...
    'nx',       nx, ...
    'nA',       nA, ...
    'neqnln',   neq, ...
    'niqnln',   niq, ...
    'dgs',      dgs, ...
    'dhs',      dhs, ...
    'Hs',       Hs    );

%% define variable and constraint bounds
options.lb = xmin;
options.ub = xmax;
options.cl = [zeros(neq, 1);  -Inf(niq, 1); l];
options.cu = [zeros(neq, 1); zeros(niq, 1); u];

%% assign function handles
funcs.objective         = @objective;
funcs.gradient          = @gradient;
funcs.constraints       = @constraints;
funcs.jacobian          = @jacobian;
funcs.hessian           = @hessian;
funcs.jacobianstructure = @(d) Js;
funcs.hessianstructure  = @(d) tril(Hs);

%%-----  run solver  -----
%% run the optimization
if have_feature('ipopt_auxdata')
    [x, info] = ipopt_auxdata(x0,funcs,options);
else
    [x, info] = ipopt(x0,funcs,options);
end

if info.status == 0 || info.status == 1
    eflag = 1;
else
    eflag = 0;
end
output = struct('status', info.status);
if isfield(info, 'iter')
    output.iterations = info.iter;
else
    output.iterations = [];
end
if isfield(info, 'cpu')
    output.cpu = info.cpu;
end
if isfield(info, 'eval')
    output.eval = info.eval;
end
f = f_fcn(x);

%% check for empty results (in case optimization failed)
if isempty(info.lambda)
    lam = NaN(nA, 1);
else
    lam = info.lambda;
end
if isempty(info.zl) && ~isempty(xmin)
    zl = NaN(nx, 1);
else
    zl = info.zl;
end
if isempty(info.zu) && ~isempty(xmax)
    zu = NaN(nx, 1);
else
    zu = info.zu;
end

%% extract shadow prices for equality var bounds converted to eq constraints
%% (since IPOPT does not return shadow prices on variables that it eliminates)
if nk
    offset = neq + niq + nA - nk;
    lam_tmp = lam(offset+(1:nk));
    kl = find(lam_tmp < 0);         %% lower bound binding
    ku = find(lam_tmp > 0);         %% upper bound binding
    zl(kk(kl)) = -lam_tmp(kl);
    zu(kk(ku)) =  lam_tmp(ku);
    lam(offset+(1:nk)) = [];        %% remove these shadow prices
    nA = nA - nk;                   %% reduce dimension accordingly
end

%% extract multipliers for linear constraints
lam_lin = lam(neq+niq+(1:nA));      %% lambda for linear constraints
kl = find(lam_lin < 0);             %% lower bound binding
ku = find(lam_lin > 0);             %% upper bound binding
mu_l = zeros(nA, 1);
mu_l(kl) = -lam_lin(kl);
mu_u = zeros(nA, 1);
mu_u(ku) = lam_lin(ku);

lambda = struct( ...
    'lower', zl, ...
    'upper', zu, ...
    'eqnonlin', lam(1:neq), ...
    'ineqnonlin', lam(neq+(1:niq)), ...
    'mu_l', mu_l, ...
    'mu_u', mu_u );


%-----  callback functions  -----
function f = objective(x, d)
f = d.f_fcn(x);

function df = gradient(x, d)
[f, df] = d.f_fcn(x);

function c = constraints(x, d)
if isempty(d.gh_fcn)
    c = d.A*x;
else
    [h, g] = d.gh_fcn(x);
    c = [g; h; d.A*x];
end

function J = jacobian(x, d)
if isempty(d.gh_fcn)
    J = d.A;
else
    [h, g, dh, dg] = d.gh_fcn(x);
    %% add sparse structure (with tiny values) to current matrices to
    %% ensure that sparsity structure matches that supplied
    J = [(dg + d.dgs)'; (dh + d.dhs)'; d.A];
end

function H = hessian(x, sigma, lambda, d)
if isempty(d.hess_fcn)
    [f_, df_, d2f] = d.f_fcn(x);    %% cost
    if ~issparse(d2f), d2f = sparse(d2f); end   %% convert to sparse if necessary
    H = tril(d2f * sigma);
else
    lam.eqnonlin   = lambda(1:d.neqnln);
    lam.ineqnonlin = lambda(d.neqnln+(1:d.niqnln));
    H = d.hess_fcn(x, lam, sigma);
    if ~issparse(H), H = sparse(H); end     %% convert to sparse if necessary
    %% add sparse structure (with tiny values) to current matrices to
    %% ensure that sparsity structure matches that supplied
    H = tril(H + d.Hs);
end
