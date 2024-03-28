function [x, f, eflag, varargout] = nleqs_fsolve(fcn, x0, opt)
% nleqs_fsolve-  Nonlinear Equation Solver based on Newton's method.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, JAC] = NLEQS_FSOLVE(FCN, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT, JAC] = NLEQS_FSOLVE(PROBLEM)
%   A wrapper function providing a standardized interface for using
%   FSOLVE to solve the nonlinear equation f(x) = 0, beginning from a
%   starting point x0.
%
%   Inputs:
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
%           max_it (0) - maximum number of iterations
%                       (0 means use solver's own default)
%           tol (0) - tolerance on f(x)
%                       (0 means use solver's own default)
%           fsolve_opt - options struct for FSOLVE, values to be passed
%                   directly to OPTIMSET or OPTIMOPTIONS, values in verbose
%                   max_it, tol override these options
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
%           alg - algorithm code of solver used ('FSOLVE')
%           iterations - number of iterations performed
%           message - exit message
%       JAC : final Jacobian matrix, J(x)
%
%   Note the calling syntax is almost identical to that of FSOLVE from
%   MathWorks' Optimization Toolbox. The function for evaluating the
%   nonlinear function and Jacobian is identical.
%
%   Calling syntax options:
%       [x, f, exitflag, output, jac] = nleqs_fsolve(fcn, x0);
%       [x, f, exitflag, output, jac] = nleqs_fsolve(fcn, x0, opt);
%       x = nleqs_fsolve(problem);
%               where problem is a struct with fields: fcn, x0, opt
%               and all fields except 'fcn' and 'x0' are optional
%       x = nleqs_fsolve(...);
%       [x, f] = nleqs_fsolve(...);
%       [x, f, exitflag] = nleqs_fsolve(...);
%       [x, f, exitflag, output] = nleqs_fsolve(...);
%       [x, f, exitflag, output, jac] = nleqs_fsolve(...);
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
%       [x, f, exitflag, output, jac] = nleqs_fsolve(problem);
%
% See also nleqs_master, fsolve.

%   MP-Opt-Model
%   Copyright (c) 2020-2024, Power Systems Engineering Research Center (PSERC)
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
        opt = struct();
    end
end

%% set default options
if isfield(opt, 'verbose') && ~isempty(opt.verbose)
    verbose = opt.verbose;
else
    verbose = 0;
end
if isfield(opt, 'fsolve_opt')
    fsolve_opt = opt.fsolve_opt;
else
    fsolve_opt = struct();
end
if verbose > 1
    fsolve_opt.Display = 'iter';
elseif verbose > 0
    fsolve_opt.Display = 'final';
else
    fsolve_opt.Display = 'off';
end
if isfield(opt, 'max_it') && opt.max_it     %% not empty or zero
    fsolve_opt.MaxIter = opt.max_it;
end
if isfield(opt, 'tol') && opt.tol           %% not empty or zero
    fsolve_opt.TolFun = opt.tol;
end
fsolve_opt.Jacobian = 'on';
fsolve_opt.JacobMult = [];

%% call the solver
[x, f, eflag, varargout{1:nargout-3}] = fsolve(fcn, x0, fsolve_opt);

if nargout > 3
    varargout{1}.alg = 'FSOLVE';
    varargout{1}.exitflag = eflag;
end
if eflag > 0
    eflag = 1;
end
