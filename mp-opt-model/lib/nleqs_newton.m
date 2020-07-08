function [x, F, eflag, output, J] = nleqs_newton(fcn, x0, opt)
%NLEQS_NEWTON  Nonlinear Equation Solver based on Newton's method.
%   [X, F, EXITFLAG, OUTPUT, JAC] = NLEQS_NEWTON(FCN, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT, JAC] = NLEQS_NEWTON(PROBLEM)
%   A function providing a standardized interface for using Newton's
%   method to solve the nonlinear equation F(X) = 0, beginning from a
%   starting point X0.
%
%   Inputs:
%       FCN : handle to function that evaluates the function F(X) to
%           be solved and its Jacobian, J(X). Calling syntax for this
%           function is:
%               [F, J] = FCN(X)
%           If F is M x 1 and and X is N x 1, then J is the M x N matrix
%           of partial derivatives of F w.r.t. X.
%       X0 : starting value of vector X
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           max_it (10) - maximum number of iterations for Newton's method
%           tol (1e-8) - tolerance on Inf-norm of F(X)
%           newton_opt - options struct for Newton's method, with field:
%               lin_solver ('') - linear solver passed to MPLINSOLVE to solve
%                        Newton update step
%                   [ ''      - default, '\' for small probs, 'LU3' larger  ]
%                   [ '\'     - built-in backslash operator                 ]
%                   [ 'LU'    - explicit default LU decomp & back substitutn]
%                   [ 'LU3'   - 3 output arg form of LU, Gilbert-Peierls    ]
%                   [           algorithm with approximate minimum degree   ]
%                   [           (AMD) reordering                            ]
%                   [ 'LU4'   - 4 output arg form of LU, UMFPACK solver     ]
%                   [           (same as 'LU')                              ]
%                   [ 'LU5'   - 5 output arg form of LU, UMFPACK solver,    ]
%                   [           w/row scaling                               ]
%                   [   (see MPLINSOLVE for complete list of all options)   ]
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: fcn, x0, opt
%
%   Outputs (all optional, except X):
%       X : solution vector
%       F : final function value
%       EXITFLAG : exit flag
%           1 = converged
%           0 or negative values = solver specific failure codes
%       OUTPUT : output struct with the following fields:
%           alg - algorithm code of solver used ('NEWTON')
%           iterations - number of iterations performed
%           hist - struct array with trajectories of the following:
%                   normF
%           message - exit message
%       JAC : final Jacobian matrix
%
%   Note the calling syntax is almost identical to that of FSOLVE from
%   MathWorks' Optimization Toolbox. The function for evaluating the
%   nonlinear function and Jacobian is identical.
%
%   Calling syntax options:
%       [x, f, exitflag, output, jac] = nleqs_newton(fcn, x0);
%       [x, f, exitflag, output, jac] = nleqs_newton(fcn, x0, opt);
%       x = nleqs_newton(problem);
%               where problem is a struct with fields: fcn, x0, opt
%               and all fields except 'fcn' and 'x0' are optional
%       x = nleqs_newton(...);
%       [x, f] = nleqs_newton(...);
%       [x, f, exitflag] = nleqs_newton(...);
%       [x, f, exitflag, output] = nleqs_newton(...);
%       [x, f, exitflag, output, jac] = nleqs_newton(...);
%
%   Example: (problem from https://www.chilimath.com/lessons/advanced-algebra/systems-non-linear-equations/)
%       function [f, J] = f1(x)
%       f = [  x(1)   + x(2) - 1;
%             -x(1)^2 + x(2) + 5    ];
%       if nargout > 1
%           J = [1 1; -2*x(1) 1];
%       end
%
%       problem = struct( ...
%           'fcn',    @(x)f1(x), ...
%           'x0',       [0; 0], ...
%           'opt',      struct('verbose', 2) ...
%       );
%       [x, f, exitflag, output, jac] = nleqs_newton(problem);
%
%   See also NLEQS_MASTER.

%   MP-Opt-Model
%   Copyright (c) 1996-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%%----- input argument handling  -----
%% gather inputs
if nargin == 1 && isstruct(fcn) %% problem struct
    p = fcn;
    fcn = p.fcn;
    x0 = p.x0;
    if isfield(p, 'opt'),   opt = p.opt;    else,   opt = [];   end
else                            %% individual args
    if nargin < 3
        opt = [];
    end
end
nx = size(x0, 1);           %% number of variables

%% set default options
opt0 = struct(  'verbose', 0, ...
                'max_it', 10, ...
                'tol', 1e-8 );
if isempty(opt)
    opt = opt0;
end
if isfield(opt, 'verbose') && ~isempty(opt.verbose)
    verbose = opt.verbose;
else
    verbose = opt0.verbose;
end
if isfield(opt, 'max_it') && opt.max_it     %% not empty or zero
    max_it = opt.max_it;
else
    max_it = opt0.max_it;
end
if isfield(opt, 'tol') && opt.tol           %% not empty or zero
    tol = opt.tol;
else
    tol = opt0.tol;
end
if isfield(opt, 'newton_opt') && isfield(opt.newton_opt, 'lin_solver')
    lin_solver = opt.newton_opt.lin_solver;
else
    lin_solver = '';
end

%% initialize
eflag = 0;
i = 0;
x = x0;
hist(max_it+1) = struct('normF', 0);

%% evaluate F(x0)
[F, J] = fcn(x);

%% check tolerance
normF = norm(F, inf);
if verbose > 1
    fprintf('\n it     max residual');
    fprintf('\n----  ----------------');
    fprintf('\n%3d     %10.3e', i, normF);
end
if normF < tol
    eflag = 1;
    msg = sprintf('Newton''s method converged in %d iterations.', i);
    if verbose > 1
        fprintf('\nConverged!\n');
    end
end

%% attempt to pick fastest linear solver, if not specified
if isempty(lin_solver)
    nf = length(F);
    if nf <= 10 || have_fcn('octave')
        lin_solver = '\';       %% default \ operator
    else    %% MATLAB and nf > 10
        lin_solver = 'LU3';     %% LU decomp with 3 output args, AMD ordering
    end
end

%% save history
hist(i+1).normF = normF;

%% do Newton iterations
while (~eflag && i < max_it)
    %% update iteration counter
    i = i + 1;

    %% compute update step
    dx = mplinsolve(J, -F, lin_solver);

    %% update voltage
    x = x + dx;

    %% evalute F(x) and J(x)
    [F, J] = fcn(x);

    %% check for convergence
    normF = norm(F, inf);
    if verbose > 1
        fprintf('\n%3d     %10.3e', i, normF);
    end

    %% save history
    hist(i+1).normF = normF;

    if normF < tol
        eflag = 1;
        msg = sprintf('Newton''s method converged in %d iterations.', i);
    end
end
if eflag ~= 1
    msg = sprintf('Newton''s method did not converge in %d iterations.', i);
end
if verbose
    fprintf('\n%s\n', msg);
end
if nargout > 3
    output = struct('alg', 'NEWTON', ...
                    'iterations', i, ...
                    'hist', hist(1:i+1), ...
                    'message', msg  );
end
