function [x, f, eflag, output, lambda] = qps_mips(H, c, A, l, u, xmin, xmax, x0, opt)
%QPS_MIPS  Quadratic Program Solver based on MIPS
%   Uses the Matlab Interior Point Solver (MIPS) to solve the following
%   QP (quadratic programming) problem:
%
%       min 1/2 x'*H*x + c'*x
%        x
%
%   subject to
%
%       l <= A*x <= u       (linear constraints)
%       xmin <= x <= xmax   (variable bounds)
%
%   [x, f, exitflag, output, lambda] = ...
%       qps_mips(H, c, A, l, u, xmin, xmax, x0, opt)
%
%   x = qps_mips(H, c, A, l, u)
%   x = qps_mips(H, c, A, l, u, xmin, xmax)
%   x = qps_mips(H, c, A, l, u, xmin, xmax, x0)
%   x = qps_mips(H, c, A, l, u, xmin, xmax, x0, opt)
%   x = qps_mips(problem), where problem is a struct with fields:
%                       H, c, A, l, u, xmin, xmax, x0, opt
%                       all fields except 'H', 'c', 'A' and 'l' are optional
%   x = qps_mips(...)
%   [x, f] = qps_mips(...)
%   [x, f, exitflag] = qps_mips(...)
%   [x, f, exitflag, output] = qps_mips(...)
%   [x, f, exitflag, output, lambda] = qps_mips(...)
%
%   Inputs:
%       H : matrix (possibly sparse) of quadratic cost coefficients
%       c : vector of linear cost coefficients
%       A, l, u : define the optional linear constraints. Default
%           values for the elements of l and u are -Inf and Inf,
%           respectively.
%       xmin, xmax : optional lower and upper bounds on the
%           x variables, defaults are -Inf and Inf, respectively.
%       x0 : optional starting value of optimization vector x
%       opt : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           feastol (1e-6) - termination tolerance for feasibility
%               condition
%           gradtol (1e-6) - termination tolerance for gradient
%               condition
%           comptol (1e-6) - termination tolerance for complementarity
%               condition
%           costtol (1e-6) - termination tolerance for cost condition
%           max_it (150) - maximum number of iterations
%           step_control (0) - set to 1 to enable step-size control
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
%           described above: H, c, A, l, u, xmin, xmax, x0, opt
%
%   Outputs:
%       x : solution vector
%       f : final objective function value
%       exitflag : exit flag,
%           1 = first order optimality conditions satisfied
%           0 = maximum number of iterations reached
%           -1 = numerically failed
%       output : structure with fields:
%           iterations - number of iterations performed
%           hist - struct array with trajectories of the following:
%                   feascond, gradcond, compcond, costcond, gamma,
%                   stepsize, obj, alphap, alphad
%           message - exit message
%       lambda : struct containing the Langrange and Kuhn-Tucker
%           multipliers on the constraints, with fields:
%           mu_l - lower bound on linear constraints
%           mu_u - upper bound on linear constraints
%           lower - lower bound on optimization variables
%           upper - upper bound on optimization variables
%
%   Note the calling syntax is almost identical to that of 'quadprog'
%   from MathWorks' Optimization Toolbox. The main difference is that
%   the linear constraints are specified with A, l, u instead of
%   A, b, Aeq, beq.
%
%   See also MIPS

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2010 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

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
