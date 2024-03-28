function varargout = nleqs_newton(varargin)
% nleqs_newton - Nonlinear Equation Solver based on Newton's method.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, JAC] = NLEQS_NEWTON(FCN, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT, JAC] = NLEQS_NEWTON(PROBLEM)
%   A function providing a standardized interface for using Newton's
%   method to solve the nonlinear equation f(x) = 0, beginning from a
%   starting point x0.
%
%   Calls NLEQS_CORE with a Newton update function.
%
%   Inputs:
%       FCN : handle to function that evaluates the function f(x) to
%           be solved and its Jacobian, J(x). Calling syntax for this
%           function is:
%               f = FCN(x)
%               [f, J] = FCN(x)
%           If f and x are n x 1, then J is the n x n matrix of partial
%           derivatives of f (rows) w.r.t. x (cols).
%       X0 : starting value, x0, of vector x
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           max_it (10) - maximum number of iterations for Newton's method
%           tol (1e-8) - tolerance on Inf-norm of f(x)
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
%       X : solution vector x
%       F : final function value, f(x)
%       EXITFLAG : exit flag
%           1 = converged
%           0 or negative values = solver specific failure codes
%       OUTPUT : output struct with the following fields:
%           alg - algorithm code of solver used ('NEWTON')
%           iterations - number of iterations performed
%           hist - struct array with trajectories of the following:
%                   normf
%           message - exit message
%       JAC : final Jacobian matrix, J(x)
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
% See also nleqs_master, nleqs_core.

%   MP-Opt-Model
%   Copyright (c) 1996-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%%----- input argument handling  -----
%% gather inputs
if nargin == 1 && isstruct(varargin{1}) %% problem struct
    p = varargin{1};
    fcn = p.fcn;
    x0 = p.x0;
    if isfield(p, 'opt'),   opt = p.opt;    else,   opt = [];   end
else                                    %% individual args
    fcn = varargin{1};
    x0  = varargin{2};
    if nargin < 3
        opt = [];
    else
        opt = varargin{3};
    end
end

%% set default options
if isempty(opt)
    opt = struct();
end
if isfield(opt, 'newton_opt') && isfield(opt.newton_opt, 'lin_solver')
    lin_solver = opt.newton_opt.lin_solver;
else
    lin_solver = '';
end

%% attempt to pick fastest linear solver, if not specified
if isempty(lin_solver)
    nx = size(x0, 1);           %% number of variables
    if nx <= 10 || have_feature('octave')
        lin_solver = '\';       %% default \ operator
    else    %% MATLAB and nx > 10
        lin_solver = 'LU3';     %% LU decomp with 3 output args, AMD ordering
    end
end

sp = struct( ...
    'alg',              'NEWTON', ...
    'name',             'Newton''s', ...
    'default_max_it',   10, ...
    'need_jac',         1, ...
    'update_fcn',       @(x, f, J)newton_update_fcn(x, f, J, lin_solver)  );

[varargout{1:nargout}] = nleqs_core(sp, fcn, x0, opt);
% opt.alg = 'CORE';
% opt.core_sp = sp;
% [varargout{1:nargout}] = nleqs_master(fcn, x0, opt);


function x = newton_update_fcn(x, f, J, lin_solver)
dx = mplinsolve(J, -f, lin_solver);     %% compute update step
x = x + dx;                             %% update x
