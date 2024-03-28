function pne_ex1
% pne_ex1 - Example of parameterized nonlinear equation solution tracing.
%
% Example of tracing the solution of the following set of parameterized
% nonlinear equations, first using opt_model and opt_model.solve, then
% directly using pnes_master.
%
% .. math::
%
%   \rvec{f}(\x, \lambda) = \left[\begin{array}{c}
%           x_1 + x_2 + 6 \lambda - 1 \\
%           -x_1^2 + x_2 + 5
%       \end{array}\right] = 0

%% variable initial values
y0 = [-1; 0];
lam0 = 0;

%% solver options
opt = struct( ...
    'verbose', 2, ...
    'stop_at', 'FULL', ...
    'step', 0.6, ...
    'plot', struct( ...
        'level', 2, ...
        'idx', 1:2, ...
        'title2', 'PNE Continuation Example', ...
        'legend', 'x_%d' ));

%%-----  METHOD 1  -----
fprintf('\n-----  METHOD 1 -----\n');
%% build model
om = opt_model;
om.add_var('y', 2, y0);
om.add_var('lam', 1, lam0);
om.add_nln_constraint('f', 2, 1, @f1p, [], {'y', 'lam'});

%% solve model
[x, f, exitflag, output, jac] = om.solve(opt);

%% print results
fprintf('\nx = \n');
fprintf('%4g\n', x);
fprintf('\nf = \n');
fprintf('%13g\n', f);
fprintf('\njac =\n');
fprintf('%4g%4g%4g\n', full(jac'));

%%-----  METHOD 2  -----
fprintf('\n-----  METHOD 2 -----\n');
%% assemble model parameters manually
x0 = [y0; lam0];

%% solve model
[x, f, exitflag, output, jac] = pnes_master(@f1p, x0, opt);

%% print results
fprintf('\nexitflag = %d\n', exitflag);
fprintf('output.max_lam = %g\n', output.max_lam);
fprintf('\nx = \n');
fprintf('%4g\n', x);
fprintf('\nf = \n');
fprintf('%13g\n', f);
fprintf('\njac =\n');
fprintf('%4g%4g%4g\n', jac');

%% parameterized 2-d problem
%% based on https://www.chilimath.com/lessons/advanced-algebra/systems-non-linear-equations/
function [f, J] = f1p(x)
if iscell(x)    %% METHOD 1
    [y, lam] = deal(x{:});
    x = [y; lam];
end
f = [  x(1)   + x(2) + 6*x(3) - 1;
      -x(1)^2 + x(2) + 5    ];
if nargout > 1
    J = [1 1 6; -2*x(1) 1 0];
end
