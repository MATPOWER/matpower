function [x, f, eflag, output, J] = nleqs_fd_newton(fcn, x0, opt)
% nleqs_fd_newton - Nonlinear Equation Solver based on Newton's method.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, JAC] = NLEQS_FD_NEWTON(FCN, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT, JAC] = NLEQS_FD_NEWTON(PROBLEM)
%   A function providing a standardized interface for using Newton's
%   method to solve the nonlinear equation f(x) = 0, beginning from a
%   starting point x0.
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
%           max_it (30) - maximum number of iterations for fast-decoupled
%               Newton's method
%           tol (1e-8) - tolerance on Inf-norm of f(x)
%           fd_opt - options struct for Newton's method, with field:
%               jac_approx_fcn (required) - handle to function with the
%                   following calling syntax:
%                       JJ = jac_approx_fcn();
%                   where JJ is a cell array of matrices whose number of
%                   rows sum to the dimension of f and number of columns
%                   sum to the dimension of x. JJ{k} is a matrix
%                   approximating the Jacobian of block k of f and x.
%               labels (optional) - cell array of same length as JJ returned
%                   by jac_approx_fcn(), with short string labels for
%                   the decoupled blocks
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
%           alg - algorithm code of solver used ('FD')
%           iterations - number of iterations performed
%           hist - struct array with trajectories of the following:
%                   normf
%           message - exit message
%       JAC : approximate Jacobian matrix
%
%   Note the calling syntax is almost identical to that of FSOLVE from
%   MathWorks' Optimization Toolbox. The function for evaluating the
%   nonlinear function is identical.
%
%   Calling syntax options:
%       [x, f, exitflag, output, jac] = nleqs_fd_newton(fcn, x0);
%       [x, f, exitflag, output, jac] = nleqs_fd_newton(fcn, x0, opt);
%       x = nleqs_fd_newton(problem);
%               where problem is a struct with fields: fcn, x0, opt
%               and all fields except 'fcn' and 'x0' are optional
%       x = nleqs_fd_newton(...);
%       [x, f] = nleqs_fd_newton(...);
%       [x, f, exitflag] = nleqs_fd_newton(...);
%       [x, f, exitflag, output] = nleqs_fd_newton(...);
%       [x, f, exitflag, output, jac] = nleqs_fd_newton(...);
%
%   Example: (problem from Christi Patton Luks, https://www.youtube.com/watch?v=pJG4yhtgerg)
%       function f = f2(x)
%       f = [  x(1)^2 +   x(1)*x(2)   - 10;
%              x(2)   + 3*x(1)*x(2)^2 - 57  ];
%
%       function JJ = jac_approx_fcn2()
%       J = [7 2; 27 37];
%       JJ = {J(1,1), J(2,2)};
%
%       problem = struct( ...
%           'fcn', @(x)f2(x), ...
%           'x0',  [0; 0], ...
%           'opt', struct( ...
%               'verbose', 2, ...
%               'fd_opt', struct( ...
%                   'jac_approx_fcn', @jac_approx_fcn2 )));
%       [x, f, exitflag, output, jac] = nleqs_fd_newton(problem);
%
% See also nleqs_master.

%   MP-Opt-Model
%   Copyright (c) 1996-2024, Power Systems Engineering Research Center (PSERC)
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
                'max_it', 30, ...
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
if isfield(opt, 'fd_opt') && isfield(opt.fd_opt, 'jac_approx_fcn')
    jac_approx_fcn = opt.fd_opt.jac_approx_fcn;
    if isfield(opt.fd_opt, 'labels')
        labels = opt.fd_opt.labels;
    else
        labels = {};
    end
else
    error('nleqs_fd_newton: required ''fd_opt.jac_approx_fcn'' option missing');
end

%% initialize
lu_vec = have_feature('lu_vec');
eflag = 0;
i = 0;
x = x0;
hist(max_it+1) = struct('normf', 0);

%% get Jacobian approximation matrices
JJ = jac_approx_fcn();

%% get indices for each block
nj = length(JJ);
normf = zeros(nj, 1);   %% initialize
b = 0;
i1 = zeros(1, nj); iN = i1;
for j = 1:nj
    [m, n] = size(JJ{j});
    if m ~= n
        error('nleqs_fd_newton: Jacobian approximation for block %d is not square (%d x %d)', j, m, n);
    end
    i1(j) = b+1;
    iN(j) = b+m;
    b = iN(j);
end

%% block labels
if isempty(labels)  %% create default labels 'A', 'B', ...
    for j = 1:nj
        labels{j} = sprintf('%s', char('A'+j-1));
    end
elseif length(labels) ~= nj
    error('nleqs_fd_newton: size of labels must be consistent with jac_approx_fcn');
end

%% evaluate f(x0)
f = fcn(x);

%% check tolerance
for j = 1:nj
    normf(j) = norm(f(i1(j):iN(j)), inf);
end
if verbose > 1
    fprintf('\n iteration ');
    for j = 1:nj, fprintf('   max residual ');          end
    fprintf('\nblock    # ');
    for j = 1:nj
        n1 = length(labels{j});
        n2 = fix((14 - 3 - n1) / 2)+2;
        n3 = 13-n1-n2;
        lb = sprintf('%sf[%s]%s', repmat(' ', 1, n2), labels{j}, repmat(' ', 1, n3));
        fprintf('%16s', lb);
    end
    fprintf('\n------ ----');
    for j = 1:nj, fprintf('  --------------');          end
    fprintf('\n  -    %3d', i);
    for j = 1:nj, fprintf('      %10.3e', normf(j));    end
end
if max(normf) < tol
    eflag = 1;
    msg = sprintf('Fast-decoupled Newton''s method converged in %d iterations.', i);
    if verbose > 1
        fprintf('\nConverged!\n');
    end
end

%% factor Jacobian approximation matrices
if lu_vec
    for j = nj:-1:1     %% backwards to init cell arrays to full-size at start
        if size(JJ{j}, 1) > 1
            [L{j}, U{j}, p{j}, q] = lu(JJ{j}, 'vector');
            [junk, iq{j}] = sort(q);
        else
            [L{j}, U{j}, p{j}] = lu(JJ{j}, 'vector');
            iq{j} = 1;
        end
    end
else
    for j = nj:-1:1     %% backwards to init cell arrays to full-size at start
        [L{j}, U{j}, P{j}] = lu(J{j});
    end
end

%% save history
hist(i+1).normf = normf;

%% do Newton iterations
while (~eflag && i < max_it)
    %% update iteration counter
    i = i + 1;

    %% do decoupled iterations
    for j = 1:nj
        ff = f(i1(j):iN(j));
        if lu_vec
            dx = -( U{j} \  (L{j} \ ff(p{j})) );
            dx = dx(iq{j});
        else
            dx = -( U{j} \  (L{j} \ (P{j} * ff)) );
        end

        %% update x
        x(i1(j):iN(j)) = x(i1(j):iN(j)) + dx;

        %% evalute f(x) and J(x)
        f = fcn(x);

        %% check for convergence
        for jj = 1:nj
            normf(jj) = norm(f(i1(jj):iN(jj)), inf);
        end
        if verbose > 1
            n1 = length(labels{j});
            n2 = fix((6 - n1)/2);
            n3 = 6-n1-n2;
            lb = sprintf('%s%s%s', repmat(' ', 1, n2), labels{j}, repmat(' ', 1, n3));
            fprintf('\n%s %3d', lb, i);
            for jj = 1:nj
                fprintf('      %10.3e', normf(jj));
            end
        end
        if max(normf) < tol
            eflag = 1;
            msg = sprintf('Fast-decoupled Newton''s method converged in ');
            for jj = 1:nj
                msg = sprintf('%s%d %s-', msg, i - (j<jj), labels{jj});
                if jj ~= nj
                    msg = sprintf('%s and ', msg);
                end
            end
            msg = sprintf('%siterations.', msg);
            break;
        end
    end

    %% save history
    hist(i+1).normf = normf;
end
if eflag ~= 1
    msg = sprintf('Fast-decoupled Newton''s method did not converge in %d iterations.', i);
end
if verbose
    fprintf('\n%s\n', msg);
end
if nargout > 3
    output = struct('alg', 'FD', ...
                    'iterations', i, ...
                    'hist', hist(1:i+1), ...
                    'message', msg  );
    if nargout > 4
        nx = length(x);
        if issparse(JJ{1})
            J = sparse(nx, nx);
        else
            J = zeros(nx, nx);
        end
        for j = 1:nj
            J(i1(j):iN(j), i1(j):iN(j)) = JJ{j};
        end
    end
end
