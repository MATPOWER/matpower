function qp_ex1
% qp_ex1 - Example of quadratic program (QP) optimization.
%
% Example of solving the following QP problem, first using opt_model and
% opt_model.solve, then directly using qps_master.
%
% .. math:: \min_{\x} \frac{1}{2} \trans{\x} \param{\rmat{Q}} \x
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
%   \param{\rmat{Q}} = \left[\begin{array}{cccc}
%            8 &  1 & -3 & -4 \\
%            1 &  4 & -2 & -1 \\
%           -3 & -2 &  5 &  4 \\
%           -4 & -1 &  4 & 12
%       \end{array}\right]
%
% .. math::
%
%   \param{\rvec{l}} = \left[\begin{array}{c} 4 \\ -\infty \end{array}\right],
%   \param{\rmat{A}} = \left[\begin{array}{cccc}
%            6 & 1 & 5 & -4 \\
%            4 & 9 & 0 &  0
%       \end{array}\right],
%   \param{\rvec{u}} = \left[\begin{array}{c} 4 \\ 2 \end{array}\right]
%
% .. math::
%
%   \param{\x}_\mathrm{min} = \left[\begin{array}{c} 0 \\ 0 \\ -\infty \\ -\infty \end{array}\right],
%   \param{\x}_\mathrm{max} = \left[\begin{array}{c} \infty \\ \infty \\ 0 \\ 2 \end{array}\right]

%% variable initial values
y0 = [1; 0];
z0 = [0; 1];

%% variable lower bounds
ymin = [0; 0];
zmax = [0; 2];

%% constraint data
A1 = [ 6 1 5 -4 ];  b1 = 4;
A2 = [ 4 9 ];       u2 = 2;

%% quadratic cost coefficients
Q = [ 8  1 -3 -4;
      1  4 -2 -1;
     -3 -2  5  4;
     -4 -1  4  12  ];

%% solver options
opt = struct('verbose', 2, 'alg', 'MIPS');

%%-----  METHOD 1  -----
%% build model
om = opt_model;
om.add_var('y', 2, y0, ymin);
om.add_var('z', 2, z0, [], zmax);
om.add_lin_constraint('lincon1', A1, b1, b1);
om.add_lin_constraint('lincon2', A2, [], u2, {'y'});
om.add_quad_cost('cost', Q, []);

%% solve model
[x, f, exitflag, output, lambda] = om.solve();
% [x, f, exitflag, output, lambda] = om.solve(opt)

%% print results
fprintf('\n-----  METHOD 1 -----');
fprintf('\nf = %g      exitflag = %d\n', f, exitflag);
fprintf('\n             var bound shadow prices\n');
fprintf('     x     lambda.lower  lambda.upper\n');
fprintf('%8.4f  %10.4f  %12.4f\n', [x lambda.lower lambda.upper]');
fprintf('\nconstraint shadow prices\n');
fprintf('lambda.mu_l  lambda.mu_u\n');
fprintf('%8.4f  %11.4f\n', [lambda.mu_l lambda.mu_u]');

%%-----  METHOD 2  -----
%% assemble model parameters manually
xmin = [ymin; -Inf(2,1)];
xmax = [ Inf(2,1); zmax];
x0 = [y0; z0];
A = [ A1; A2 0 0];
l = [ b1; -Inf ];
u = [ b1;  u2  ];

%% solve model
[x, f, exitflag, output, lambda] = qps_master(Q, [], A, l, u, xmin, xmax, x0);
% [x, f, exitflag, output, lambda] = qps_master(Q, [], A, l, u, xmin, [], x0, opt)

%% print results
fprintf('\n-----  METHOD 2 -----');
fprintf('\nf = %g      exitflag = %d\n', f, exitflag);
fprintf('\n             var bound shadow prices\n');
fprintf('     x     lambda.lower  lambda.upper\n');
fprintf('%8.4f  %10.4f  %12.4f\n', [x lambda.lower lambda.upper]');
fprintf('\nconstraint shadow prices\n');
fprintf('lambda.mu_l  lambda.mu_u\n');
fprintf('%8.4f  %11.4f\n', [lambda.mu_l lambda.mu_u]');
