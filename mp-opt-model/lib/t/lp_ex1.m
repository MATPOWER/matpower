function lp_ex1
% lp_ex1 - Example of linear program (LP) optimization.
%
% Example of solving the following LP problem using qps_master.
%
% .. math:: \min_{\x} \trans{\param{\rvec{c}}} \x
%
% subject to
%
% .. math:: \param{\rvec{l}} \le \param{\rmat{A}} \x \le \param{\rvec{u}}
% .. math:: \param{\x}_\mathrm{min} \le \x \le \param{\x}_\mathrm{max}
%
% where
%
% .. math::
%
%   \param{\rvec{c}} = \left[\begin{array}{c} -5 \\ -4 \\ -6 \end{array}\right]
%
% .. math::
%
%   \param{\rvec{l}} = \left[\begin{array}{c} -\infty \\ -42 \\ -\infty \end{array}\right],
%   \param{\rmat{A}} = \left[\begin{array}{cccc}
%            1 & -1 & 1 \\
%           -3 & -2 & -4 \\
%            3 &  2 & 0
%       \end{array}\right],
%   \param{\rvec{u}} = \left[\begin{array}{c} 20 \\ \infty\\ 30 \end{array}\right]
%
% .. math::
%
%   \param{\x}_\mathrm{min} = \left[\begin{array}{c} 0 \\ 0 \\ 0 \end{array}\right],
%   \param{\x}_\mathrm{max} = \left[\begin{array}{c} \infty \\ \infty \\ \infty \end{array}\right]

%% variable initial values
x0 = [0; 0; 0];

%% variable lower bounds
xmin = [0; 0; 0];
xmax = [Inf; Inf; Inf];

%% constraint data
A = [1 -1  1;
    -3  -2  -4;
    3  2  0];
l = [-Inf; -42; -Inf];
u = [20; Inf; 30];

%% linear cost coefficients
c = [-5; -4; -6];

%% solver options
opt = struct('verbose', 2, 'alg', 'MIPS');

%% solve model
[x, f, exitflag, output, lambda] = qps_master([], c, A, l, u, xmin, xmax, x0, opt);

%% print results
fprintf('\n-----  ACTUAL RESULTS  -----\n');
fprintf('\nf = %g      exitflag = %d\n', f, exitflag);
fprintf('\n             var bound shadow prices\n');
fprintf('     x     lambda.lower  lambda.upper\n');
fprintf('%8.4f  %10.4f  %12.4f\n', [x lambda.lower lambda.upper]');
fprintf('\nconstraint shadow prices\n');
fprintf('lambda.mu_l  lambda.mu_u\n');
fprintf('%8.4f  %11.4f\n', [lambda.mu_l lambda.mu_u]');

% expected results
expected_x = [0; 15; 3];
expected_f = -78;
expected_lam.mu_l = [0;1.5;0];
expected_lam.mu_u = [0;0;0.5];
expected_lam.lower = [1;0;0];
expected_lam.upper = zeros(size(x));

%% print expected results
fprintf('\n-----  EXPECTED RESULTS  -----\n');
fprintf('\nf = %g\n', expected_f);
fprintf('\n             var bound shadow prices\n');
fprintf('     x     lambda.lower  lambda.upper\n');
fprintf('%8.4f  %10.4f  %12.4f\n', [expected_x expected_lam.lower expected_lam.upper]');
fprintf('\nconstraint shadow prices\n');
fprintf('lambda.mu_l  lambda.mu_u\n');
fprintf('%8.4f  %11.4f\n', [expected_lam.mu_l expected_lam.mu_u]');
