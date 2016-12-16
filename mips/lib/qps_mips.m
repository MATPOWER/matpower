function [x, f, eflag, output, lambda] = qps_mips(H, c, A, l, u, xmin, xmax, x0, opt)
%QPS_MIPS  Quadratic Program Solver based on MIPS.
%   [X, F, EXITFLAG, OUTPUT, LAMBDA] = ...
%       QPS_MIPS(H, C, A, L, U, XMIN, XMAX, X0, OPT)
%   Uses the MATPOWER Interior Point Solver (MIPS) to solve the following
%   QP (quadratic programming) problem:
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
%               function for improved conditioning.
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
%           described above: H, c, A, l, u, xmin, xmax, x0, opt
%
%   Outputs:
%       X : solution vector
%       F : final objective function value
%       EXITFLAG : exit flag
%           1 = first order optimality conditions satisfied
%           0 = maximum number of iterations reached
%           -1 = numerically failed
%       OUTPUT : output struct with the following fields:
%           iterations - number of iterations performed
%           hist - struct array with trajectories of the following:
%                   feascond, gradcond, compcond, costcond, gamma,
%                   stepsize, obj, alphap, alphad
%           message - exit message
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
%           qps_mips(H, c, A, l, u, xmin, xmax, x0, opt)
%
%       x = qps_mips(H, c, A, l, u)
%       x = qps_mips(H, c, A, l, u, xmin, xmax)
%       x = qps_mips(H, c, A, l, u, xmin, xmax, x0)
%       x = qps_mips(H, c, A, l, u, xmin, xmax, x0, opt)
%       x = qps_mips(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'c', 'A' and 'l' or 'u' are optional
%       x = qps_mips(...)
%       [x, f] = qps_mips(...)
%       [x, f, exitflag] = qps_mips(...)
%       [x, f, exitflag, output] = qps_mips(...)
%       [x, f, exitflag, output, lambda] = qps_mips(...)
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
%       [x, f, s, out, lambda] = qps_mips(H, c, A, l, u, xmin, [], x0, opt);
%
%   See also MIPS.

%   MIPS
%   Copyright (c) 2010-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MIPS.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mips for more info.

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
        error('qps_mips: LP problem must include constraints or variable bounds');
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
else
    nx = size(p.H, 1);
end
if ~isfield(p, 'c') || isempty(p.c)
    p.c = zeros(nx, 1);
end
if ~isfield(p, 'x0') || isempty(p.x0)
    p.x0 = zeros(nx, 1);
end

%%-----  run optimization  -----
p.f_fcn = @(x)qp_f(x, p.H, p.c);
[x, f, eflag, output, lambda] = mips(p);

%%-----  objective function  -----
function [f, df, d2f] = qp_f(x, H, c)
f = 0.5 * x' * H * x + c' * x;
if nargout > 1
    df = H * x + c;
    if nargout > 2
        d2f = H;
    end
end
