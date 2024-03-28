function varargout = nleqs_gauss_seidel(varargin)
% nleqs_gauss_seidel - Nonlinear Equation Solver based on Gauss-Seidel method.
% ::
%
%   [X, F, EXITFLAG, OUTPUT] = NLEQS_GAUSS_SEIDEL(FCN, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT] = NLEQS_GAUSS_SEIDEL(PROBLEM)
%   A function providing a standardized interface for using Gauss-Seidel
%   method to solve the nonlinear equation f(x) = 0, beginning from a
%   starting point x0.
%
%   Calls NLEQS_CORE with a user-provided update function implementing
%   the Gauss-Seidel update.
%
%   Inputs:
%       FCN : handle to function that evaluates the function f(x) to
%           be solved. Calling syntax for this function is:
%               f = FCN(x)
%       X0 : starting value, x0, of vector x
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           max_it (1000) - maximum number of iterations for Gauss-Seidel method
%           tol (1e-8) - tolerance on Inf-norm of f(x)
%           gs_opt - options struct for Gauss-Seidel method, with field:
%               x_update_fcn (required) - handle to function that performs
%                   the Gauss-Seidel update step, with the following calling
%                   syntax:  x = x_update_fcn(x, f);
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
%           alg - algorithm code of solver used ('GS')
%           iterations - number of iterations performed
%           hist - struct array with trajectories of the following:
%                   normf
%           message - exit message
%
%   Note the calling syntax is almost identical to that of FSOLVE from
%   MathWorks' Optimization Toolbox. The function for evaluating the
%   nonlinear function is identical.
%
%   Calling syntax options:
%       [x, f, exitflag, output] = nleqs_gauss_seidel(fcn, x0);
%       [x, f, exitflag, output] = nleqs_gauss_seidel(fcn, x0, opt);
%       x = nleqs_gauss_seidel(problem);
%               where problem is a struct with fields: fcn, x0, opt
%               and all fields except 'fcn' and 'x0' are optional
%       x = nleqs_gauss_seidel(...);
%       [x, f] = nleqs_gauss_seidel(...);
%       [x, f, exitflag] = nleqs_gauss_seidel(...);
%       [x, f, exitflag, output] = nleqs_gauss_seidel(...);
%
%   Example: (problem from Christi Patton Luks, https://www.youtube.com/watch?v=pJG4yhtgerg)
%       function f = f2(x)
%       f = [  x(1)^2 +   x(1)*x(2)   - 10;
%              x(2)   + 3*x(1)*x(2)^2 - 57  ];
%
%       function x = x_update_fcn2(x, f)
%       x(1) = sqrt(10 - x(1)*x(2));
%       x(2) = sqrt((57-x(2))/3/x(1));
%
%       problem = struct( ...
%           'fcn', @(x)f2(x), ...
%           'x0',  [0; 0], ...
%           'opt', struct( ...
%               'verbose', 2, ...
%               'gs_opt', struct( ...
%                   'x_update_fcn', @x_update_fcn2 )));
%       [x, f, exitflag, output] = nleqs_gauss_seidel(problem);
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
if ~isfield(opt, 'gs_opt') || ~isfield(opt.gs_opt, 'x_update_fcn')
    error('nleqs_gauss_seidel: required ''gs_opt.x_update_fcn'' option missing');
end

sp = struct( ...
    'alg',              'GS', ...
    'name',             'Gauss-Seidel', ...
    'default_max_it',   1000, ...
    'need_jac',         0, ...
    'update_fcn',       opt.gs_opt.x_update_fcn  );

[varargout{1:nargout}] = nleqs_core(sp, fcn, x0, opt);
% opt.alg = 'CORE';
% opt.core_sp = sp;
% [varargout{1:nargout}] = nleqs_master(fcn, x0, opt);
