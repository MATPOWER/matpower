function qp_ex1
% qp_ex1 - Example of quadratic program (QP) optimization.
%
% Example of solving the following QP problem, first using mp.opt_model and
% mp.opt_model.solve, then directly using qps_master.
%
% .. math:: \min_{\x} \frac{1}{2} \trans{\x} \Hh \x
%
% subject to
%
% .. math:: \l \le \AA \x \le \u
% .. math:: \param{\x}_\mathrm{min} \le \x \le \param{\x}_\mathrm{max}
%
% where
%
% .. math::
%
%   \Hh = \left[\begin{array}{cccc}
%            8 &  1 & -3 & -4 \\
%            1 &  4 & -2 & -1 \\
%           -3 & -2 &  5 &  4 \\
%           -4 & -1 &  4 & 12
%       \end{array}\right]
%
% .. math::
%
%   \l = \left[\begin{array}{c} 4 \\ -\infty \end{array}\right],
%   \AA = \left[\begin{array}{cccc}
%            6 & 1 & 5 & -4 \\
%            4 & 9 & 0 &  0
%       \end{array}\right],
%   \u = \left[\begin{array}{c} 4 \\ 2 \end{array}\right]
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
H = [ 8  1 -3 -4;
      1  4 -2 -1;
     -3 -2  5  4;
     -4 -1  4  12  ];

%% solver options
opt = struct('verbose', 2, 'alg', 'MIPS');

%%-----  METHOD 1  -----
%% build model
mm = mp.opt_model;
mm.var.add('y', 2, y0, ymin);
mm.var.add('z', 2, z0, [], zmax);
mm.lin.add(mm.var, 'lincon1', A1, b1, b1);
mm.lin.add(mm.var, 'lincon2', A2, [], u2, {'y'});
mm.qdc.add(mm.var, 'cost', H, []);

%% solve model
[x, f, exitflag, output, lambda] = mm.solve();
% [x, f, exitflag, output, lambda] = mm.solve(opt)

%% print results
fprintf('\n-----  METHOD 1 -----');
fprintf('\nf = %g      exitflag = %d\n', f, exitflag);
mm.display_soln();
% fprintf('\n             var bound shadow prices\n');
% fprintf('     x     lambda.lower  lambda.upper\n');
% fprintf('%8.4f  %10.4f  %12.4f\n', [x lambda.lower lambda.upper]');
% fprintf('\nconstraint shadow prices\n');
% fprintf('lambda.mu_l  lambda.mu_u\n');
% fprintf('%8.4f  %11.4f\n', [lambda.mu_l lambda.mu_u]');

%%-----  METHOD 2  -----
%% assemble model parameters manually
xmin = [ymin; -Inf(2,1)];
xmax = [ Inf(2,1); zmax];
x0 = [y0; z0];
A = [ A1; A2 0 0];
l = [ b1; -Inf ];
u = [ b1;  u2  ];

%% solve model
[x, f, exitflag, output, lambda] = qps_master(H, [], A, l, u, xmin, xmax, x0);
% [x, f, exitflag, output, lambda] = qps_master(H, [], A, l, u, xmin, [], x0, opt)

%% print results
fprintf('\n-----  METHOD 2 -----');
fprintf('\nf = %g      exitflag = %d\n', f, exitflag);
fprintf('\n             var bound shadow prices\n');
fprintf('     x     lambda.lower  lambda.upper\n');
fprintf('%8.4f  %10.4f  %12.4f\n', [x lambda.lower lambda.upper]');
fprintf('\nconstraint shadow prices\n');
fprintf('lambda.mu_l  lambda.mu_u\n');
fprintf('%8.4f  %11.4f\n', [lambda.mu_l lambda.mu_u]');
