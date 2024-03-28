function nleqs_master_ex1(alg)
% nleqs_master_ex1 - Example 1 of nonlinear equation solving.
% ::
%
%   nleqs_master_ex1()
%   nleqs_master_ex1('NEWTON')
%
% Example of solving the following set of nonlinear equations using
% nleqs_master.
%
% .. math::
%
%   \rvec{f}(\x) = \left[\begin{array}{c}
%           x_1 + x_2 - 1 \\
%           -x_1^2 + x_2 + 5
%       \end{array}\right] = 0

if nargin < 1
    alg = 'DEFAULT';
end
problem = struct( ...
    'fcn',  @f1, ...
    'x0',   [0; 0], ...
    'opt',  struct('verbose', 2, 'alg', alg) ...
);
[x, f, exitflag, output, jac] = nleqs_master(problem);
fprintf('\nexitflag = %d\n', exitflag);
fprintf('\nx = \n');
fprintf('   %2g\n', x);
fprintf('\nf = \n');
fprintf('   %12g\n', f);
fprintf('\njac =\n');
fprintf('   %2g  %2g\n', jac');

function [f, J] = f1(x)
%% from https://www.chilimath.com/lessons/advanced-algebra/systems-non-linear-equations/
f = [  x(1)   + x(2) - 1;
      -x(1)^2 + x(2) + 5    ];
if nargout > 1
    J = [1 1; -2*x(1) 1];
end
