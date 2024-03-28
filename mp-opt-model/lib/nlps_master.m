function [x, f, eflag, output, lambda] = nlps_master(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt)
% nlps_master - Nonlinear programming (NLP) Solver wrapper function.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       NLPS_MASTER(F_FCN, X0, A, L, U, XMIN, XMAX, GH_FCN, HESS_FCN, OPT)
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = NLPS_MASTER(PROBLEM)
%   A common wrapper function for various NLP solvers.
%   Solves the following NLP (nonlinear programming) problem:
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
%           alg ('DEFAULT') : determines which solver to use
%               'DEFAULT' : automatic, current default is MIPS
%               'MIPS'    : MIPS, MATPOWER Interior Point Solver
%                        pure MATLAB implementation of a primal-dual
%                        interior point method, if mips_opt.step_control = 1
%                        it uses MIPS-sc, a step controlled variant of MIPS
%               'FMINCON' : FMINCON, MATLAB Optimization Toolbox
%               'IPOPT'   : IPOPT, requires MEX interface to IPOPT solver
%                           https://github.com/coin-or/Ipopt
%               'KNITRO'  : Artelys Knitro, requires Artelys Knitro solver
%                           https://www.artelys.com/solvers/knitro/
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           mips_opt    - options struct for MIPS
%           fmincon_opt - options struct for FMINCON
%           ipopt_opt   - options struct for IPOPT
%           knitro_opt  - options struct for Artelys Knitro
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
%           0 or negative values = solver specific failure codes
%       OUTPUT : output struct with the following fields:
%           alg - algorithm code of solver used
%           (others) - algorithm specific fields
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
%           nlps_master(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
%
%       x = nlps_master(f_fcn, x0);
%       x = nlps_master(f_fcn, x0, A, l);
%       x = nlps_master(f_fcn, x0, A, l, u);
%       x = nlps_master(f_fcn, x0, A, l, u, xmin);
%       x = nlps_master(f_fcn, x0, A, l, u, xmin, xmax);
%       x = nlps_master(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn);
%       x = nlps_master(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn);
%       x = nlps_master(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
%       x = nlps_master(problem);
%               where problem is a struct with fields:
%                   f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt
%                   all fields except 'f_fcn' and 'x0' are optional
%       x = nlps_master(...);
%       [x, f] = nlps_master(...);
%       [x, f, exitflag] = nlps_master(...);
%       [x, f, exitflag, output] = nlps_master(...);
%       [x, f, exitflag, output, lambda] = nlps_master(...);
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
%       [x, f, exitflag, output, lambda] = nlps_master(problem);
%
% See also mips, nlps_fmincon, nlps_ipopt, nlps_knitro.

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

%% default options
if ~isempty(opt) && isfield(opt, 'alg') && ~isempty(opt.alg)
    alg = opt.alg;
else
    alg = 'DEFAULT';
end
if ~isempty(opt) && isfield(opt, 'verbose') && ~isempty(opt.verbose)
    verbose = opt.verbose;
else
    verbose = 0;
end
if strcmp(alg, 'DEFAULT')
    alg = 'MIPS';
end

%%----- call the appropriate solver  -----
switch alg
    case 'MIPS'                 %% use MIPS
        %% set up options
        if ~isempty(opt) && isfield(opt, 'mips_opt') && ~isempty(opt.mips_opt)
            mips_opt = opt.mips_opt;
        else
            mips_opt = [];
        end
        mips_opt.verbose = verbose;
        
        %% call solver
        [x, f, eflag, output, lambda] = ...
            mips(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, mips_opt);
    case 'FMINCON'              %% use fmincon
        [x, f, eflag, output, lambda] = ...
            nlps_fmincon(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
    case 'IPOPT'                %% use IPOPT
        [x, f, eflag, output, lambda] = ...
            nlps_ipopt(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
    case 'KNITRO'               %% use Artelys Knitro
        [x, f, eflag, output, lambda] = ...
            nlps_knitro(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
    otherwise
        fcn = ['nlps_' lower(alg)];
        if exist([fcn '.m']) == 2
            opt_name = [lower(alg) '_opt'];
            if isfield(opt, opt_name)
                alg_opt = opt.(opt_name);
            else
                alg_opt = struct();
            end
            alg_opt.verbose = verbose;
            [x, f, eflag, output, lambda] = ...
                feval(fcn, f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, alg_opt);
        else
            error('nlps_master: ''%s'' is not a valid algorithm code', alg);
        end
end
if ~isfield(output, 'alg') || isempty(output.alg)
    output.alg = alg;
end
