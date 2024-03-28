function nleqs_master_ex2(alg)
% nleqs_master_ex2 - Example 2 of nonlinear equation solving.
% ::
%
%   nleqs_master_ex2()
%   nleqs_master_ex2('NEWTON')
%   nleqs_master_ex2('GS')
%   nleqs_master_ex2('FD')
%
% Example of solving the following set of nonlinear equations using
% nleqs_master.
%
% .. math::
%
%   \rvec{f}(\x) = \left[\begin{array}{c}
%           x_1^2 + x_1 x_2 - 10 \\
%           x_2 + 3 x_1 x_2^2 - 57
%       \end{array}\right] = 0

if nargin < 1
    alg = 'DEFAULT';
end
x0 = [1; 2];
opt = struct( ...
    'verbose', 2, ...
    'alg', alg, ...
    'fd_opt', struct( ...
        'jac_approx_fcn', @jac_approx_fcn2, ...
        'labels', {{'P','Q'}}), ...
    'gs_opt', struct('x_update_fcn', @x_update_fcn2) );
[x, f, exitflag, output] = nleqs_master(@f2, x0, opt);
fprintf('\nexitflag = %d\n', exitflag);
fprintf('\nx = \n');
fprintf('   %2g\n', x);
fprintf('\nf = \n');
fprintf('   %12g\n', f);

function [f, J] = f2(x)
%% from Christi Patton Luks, https://www.youtube.com/watch?v=pJG4yhtgerg
f = [  x(1)^2 +   x(1)*x(2)   - 10;
       x(2)   + 3*x(1)*x(2)^2 - 57  ];
if nargout > 1
    J = [   2*x(1)+x(2)    x(1);
            3*x(2)^2       6*x(1)*x(2)+1    ];
end

function JJ = jac_approx_fcn2()
%% for use with fast-decoupled Newton's method
J = [7 2; 27 37];
JJ = {J(1,1), J(2,2)};

function x = x_update_fcn2(x, f)
%% for use with Gauss-Seidel method
x(1) = sqrt(10 - x(1)*x(2));
x(2) = sqrt((57-x(2))/3/x(1));
