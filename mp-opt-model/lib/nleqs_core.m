function [x, f, eflag, output, J] = nleqs_core(sp, fcn, x0, opt)
% nleqs_core - Core Nonlinear Equation Solver with arbitrary update function.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, JAC] = NLEQS_CORE(SP, FCN, X0, OPT)
%   A function providing the core code for a standardized interface for
%   various methods of solving the nonlinear equation f(x) = 0, beginning
%   from a starting point x0. Allows for an arbitrary update function.
%
%   Inputs:
%       SP : solver parameters, struct with the following fields (all required)
%           alg : algorithm code, e.g. 'NEWTON', 'GS', etc.
%           name : user-visible name to be followed by 'method' in output
%               'Newton''s', 'Gauss-Seidel', etc.
%           default_max_it : default value for max_it option, if not provided
%           need_jac : 1 - compute Jacobian, 0 - do not compute Jacobian
%           update_fcn : handle to function used to update x, with syntax:
%               x = update_fcn(x, f)
%               x = update_fcn(x, f, J)
%       FCN : handle to function that evaluates the function f(x) to
%           be solved and (optionally) its Jacobian, J(x). Calling syntax
%           for this function is:
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
%           max_it (SP.default_max_it) - maximum number of iterations
%           tol (1e-8) - tolerance on Inf-norm of f(x)
%
%   Outputs (all optional, except X):
%       X : solution vector x
%       F : final function value, f(x)
%       EXITFLAG : exit flag
%           1 = converged
%           0 or negative values = solver specific failure codes
%       OUTPUT : output struct with the following fields:
%           alg - algorithm code of solver used (SP.alg)
%           iterations - number of iterations performed
%           hist - struct array with trajectories of the following:
%                   normf
%                   normdx
%           message - exit message
%       JAC : final Jacobian matrix, J
%
%   Calling syntax options:
%       [x, f, exitflag, output, jac] = nleqs_core(sp, fcn, x0);
%       [x, f, exitflag, output, jac] = nleqs_core(sp, fcn, x0, opt);
%       x = nleqs_core(...);
%       [x, f] = nleqs_core(...);
%       [x, f, exitflag] = nleqs_core(...);
%       [x, f, exitflag, output] = nleqs_core(...);
%       [x, f, exitflag, output, jac] = nleqs_core(...);
%
%   Example: (problem from https://www.chilimath.com/lessons/advanced-algebra/systems-non-linear-equations/)
%       function [f, J] = f1(x)
%       f = [  x(1)   + x(2) - 1;
%             -x(1)^2 + x(2) + 5    ];
%       if nargout > 1
%           J = [1 1; -2*x(1) 1];
%       end
%
%       fcn = @(x)f1(x);
%       x0  = [0; 0];
%       opt = struct('verbose', 2);
%       sp = struct( ...
%           'alg',              'NEWTON', ...
%           'name',             'Newton''s', ...
%           'default_max_it',   10, ...
%           'need_jac',         1, ...
%           'update_fcn',       @newton_update_fcn, ...
%           )
%       [x, f, exitflag, output, jac] = nleqs_core(sp, fcn, x0, opt);
%
% See also nleqs_master, nleqs_newton, nleqs_gauss_seidel.

%   MP-Opt-Model
%   Copyright (c) 1996-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% set default options
opt0 = struct(  'verbose', 0, ...
                'max_it', sp.default_max_it, ...
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

%% initialize
eflag = 0;
i = 0;
x = x0;
hist(max_it+1) = struct('normf', 0, 'normdx', 0);

%% evaluate f(x0)
if sp.need_jac
    [f, J] = fcn(x);
else
    f = fcn(x);
end

%% check tolerance
normf = norm(f, inf);
if verbose > 1
    fprintf('\n it    max residual        max ∆x');
    fprintf('\n----  --------------  --------------');
    fprintf('\n%3d     %10.3e           -    ', i, normf);
end
if normf < tol
    eflag = 1;
    msg = sprintf('%s method converged in %d iterations.', sp.name, i);
    if verbose > 1
        fprintf('\nConverged!\n');
    end
end

%% save history
hist(i+1).normf = normf;

%% do solver iterations
while (~eflag && i < max_it)
    %% update iteration counter
    i = i + 1;

    %% save previous x
    xp = x;

    if sp.need_jac
        %% update x
        x = sp.update_fcn(x, f, J);

        %% evalute f(x) and J(x)
        [f, J] = fcn(x);
    else
        %% update x
        x = sp.update_fcn(x, f);

        %% evaluate f(x)
        f = fcn(x);
    end

    %% check for convergence
    normf  = norm(f, inf);
    normdx = norm(x-xp, inf);
    if verbose > 1
        fprintf('\n%3d     %10.3e      %10.3e', i, normf, normdx);
    end
    if normf < tol
        eflag = 1;
        msg = sprintf('%s method converged in %d iterations.', sp.name, i);
    end

    %% save history
    hist(i+1).normf  = normf;
    hist(i+1).normfx = normdx;
end
if eflag ~= 1
    msg = sprintf('%s method did not converge in %d iterations.', sp.name, i);
end
if verbose
    fprintf('\n%s\n', msg);
end
if nargout > 3
    output = struct('alg', sp.alg, ...
                    'iterations', i, ...
                    'hist', hist(1:i+1), ...
                    'message', msg  );
    if nargout > 4 && ~sp.need_jac
        J = [];
    end
end
