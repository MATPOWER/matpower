function varargout = nleqs_master(fcn, x0, opt)
% nleqs_master - Nonlinear Equation Solver wrapper function.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, JAC] = NLEQS_MASTER(FCN, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT, JAC] = NLEQS_MASTER(PROBLEM)
%   A common wrapper function for various nonlinear equation solvers.
%   Solves the nonlinear equation f(x) = 0, beginning from a starting
%   point x0.
%
%   Inputs:
%       FCN : handle to function that evaluates the function f(x) to
%           be solved and its (optionally, depending on the selected
%           solver) Jacobian, J(x). Calling syntax for this function is:
%               f = FCN(x)
%               [f, J] = FCN(x)
%           If f and x are n x 1, then J is the n x n matrix of partial
%           derivatives of f (rows) w.r.t. x (cols).
%       X0 : starting value, x0, of vector x
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           alg ('DEFAULT') : determines which solver to use
%               'DEFAULT' : automatic, current default is NEWTON
%               'NEWTON'  : standard, full-Jacobian Newton's method
%               'CORE'    : core algorithm, with arbitrary update function
%               'FD'      : fast-decoupled Newton's method
%               'FSOLVE'  : FSOLVE, MATLAB Optimization Toolbox
%               'GS'      : Gauss-Seidel method
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1 = some progress output
%               2 = verbose progress output
%           max_it (0) - maximum number of iterations
%                       (0 means use solver's own default)
%           tol (0) - termination tolerance on f(x)
%                       (0 means use solver's own default)
%           core_sp - solver parameters struct for NLEQS_CORE, required
%               when alg = 'CORE' (see NLEQS_CORE for details)
%           fd_opt - options struct for fast-decoupled Newton, NLEQS_FD_NEWTON
%           fsolve_opt - options struct for FSOLVE
%           gs_opt - options struct for Gauss-Seidel method, NLEQS_GAUSS_SEIDEL
%           newton_opt - options struct for Newton's method, NLEQS_NEWTON
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
%           alg - algorithm code of solver used
%           (others) - algorithm specific fields
%       JAC : final Jacobian matrix, J(x)
%
%   Note the calling syntax is almost identical to that of FSOLVE from
%   MathWorks' Optimization Toolbox. The function for evaluating the
%   nonlinear function and Jacobian is identical.
%
%   Calling syntax options:
%       [x, f, exitflag, output, jac] = nleqs_master(fcn, x0);
%       [x, f, exitflag, output, jac] = nleqs_master(fcn, x0, opt);
%       x = nleqs_master(problem);
%               where problem is a struct with fields: fcn, x0, opt
%               and all fields except 'fcn' and 'x0' are optional
%       x = nleqs_master(...);
%       [x, f] = nleqs_master(...);
%       [x, f, exitflag] = nleqs_master(...);
%       [x, f, exitflag, output] = nleqs_master(...);
%       [x, f, exitflag, output, jac] = nleqs_master(...);
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
%       [x, f, exitflag, output, jac] = nleqs_master(problem);
%
% See also nleqs_newton, nleqs_core, nleqs_fd_newton, nleqs_fsolve,
% nleqs_gauss_seidel, fsolve.

%   MP-Opt-Model
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
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

%% default options
if ~isempty(opt) && isfield(opt, 'alg') && ~isempty(opt.alg)
    alg = opt.alg;
else
    alg = 'DEFAULT';
end
if strcmp(alg, 'DEFAULT')
    alg = 'NEWTON';
end

%%----- call the appropriate solver  -----
switch alg
    case 'NEWTON'               %% use Newton's method solver
        nleqs_fcn = @nleqs_newton;
    case 'FD'                   %% use fast-decoupled Newton's method solver
        nleqs_fcn = @nleqs_fd_newton;
    case 'FSOLVE'               %% use fsolve
        nleqs_fcn = @nleqs_fsolve;
    case 'GS'                   %% use Gauss-Seidel solver
        nleqs_fcn = @nleqs_gauss_seidel;
    case 'CORE'                 %% use core solver
        nleqs_fcn = @(f, x, o)nleqs_core(opt.core_sp, f, x, o);
    otherwise
        error('nleqs_master: ''%s'' is not a valid algorithm code', alg);
end
[varargout{1:nargout}] = nleqs_fcn(fcn, x0, opt);
