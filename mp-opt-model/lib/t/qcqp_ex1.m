function qcqp_ex1
% qcqp_ex1 - Example of quadratically-constrained quadratic program (QCQP) optimization.
%
% Example of solving the following QCQP problem, first using opt_model and
% opt_model.solve, then directly using qcqps_master.
%
% .. math:: \min_{\x} \frac{1}{2} \trans{\x} \Hh \x + \trans{\c} \x
%
% subject to
%
% .. math::
%    :label: eq_qcqp_ex1_qc
%
%    \lq \le \frac{1}{2}\textrm{diag}\left(\trans{\diag{\{\x\}_{\times n_q}}}
%       \diag{\{\QQ_i\}_{i=1}^{n_q}} \diag{\{\x\}_{\times n_q}}\right) + \Bb \x
%       \le \uq,
% .. math:: \l \le \param{\rmat{A}} \x \le \u
% .. math:: \param{\x}_\mathrm{min} \le \x \le \param{\x}_\mathrm{max}
%
% where :eq:`eq_qcqp_ex1_qc` can also be written as
%
% .. math:: \lqi \le \frac{1}{2}\trans{\x} \QQ_i \x + \b_i \x \le \uqi, \ \ \ i = 1,2,..., n_q
%
% For the problem in this example, we have
%
% .. math::
%
%   \Hh = \left[\begin{array}{ccc}
%            0 & 0 & 0 \\
%            0 & 0 & 0 \\
%            0 & 0 & 0
%       \end{array}\right], 
%   \c = \left[\begin{array}{c}
%            0 \\
%            0 \\
%            -1
%       \end{array}\right]
%
% .. math::
%
%   \QQ_1 = \left[\begin{array}{ccc}
%            2 & 0 & 0 \\
%            0 & 2 & 0 \\
%            0 & 0 & -2
%       \end{array}\right], 
%   \QQ_2 = \left[\begin{array}{ccc}
%            2 & 0 & 0 \\
%            0 & 0 & -2 \\
%            0 & -2 & 0
%       \end{array}\right]
%
% .. math::
%
%   \lq = \left[\begin{array}{c} -\infty \\ -\infty \end{array}\right],
%   \Bb = \left[\begin{array}{ccc}
%            0 & 0 & 0 \\
%            0 & 0 & 0
%       \end{array}\right],
%   \uq = \left[\begin{array}{c} 0 \\ 0 \end{array}\right]
%
% .. math::
%
%   \l = 1, 
%   \AA = \left[\begin{array}{ccc}
%            1 & 1 & 1
%       \end{array}\right],
%   \u = 1
%
% .. math::
%
%   \param{\x}_\mathrm{min} = \left[\begin{array}{c} 0 \\ 0 \\ 0 \end{array}\right],
%   \param{\x}_\mathrm{max} = \left[\begin{array}{c} \infty \\ \infty\\ \infty
%       \end{array}\right]

%% variable initial values
y0 = [0; 0];
z0 = 1;

%% variable lower bounds
ymin = [0; 0];
zmin = 0;

%% quadratic constraint data
Q = { sparse([2 0 0; 0 2 0; 0 0 -2]); sparse([2 0 0; 0 0 -2; 0 -2 0]) };
B = zeros(2, 3);
lq = [-Inf;-Inf];
uq = [0; 0];

%% linear constraint data
A = [ 1 1 1 ];  b = 1;

%% cost data
H = [];
c = [-1;0;0];

%% solver options
opt = struct('verbose', 2, 'alg', 'MIPS');

%%-----  METHOD 1  -----
%% build model
om = opt_model().init_set_types();
om.var.add('y', 2, y0, ymin);
om.var.add('z', 1, z0, zmin);
om.qcn.add(om.var, 'quadcon', Q, B, lq, uq);
om.lin.add(om.var, 'lincon', A, b, b);
om.qdc.add(om.var, 'cost', H, c);

%% solve model
[x, f, exitflag, output, lambda] = om.solve();
% [x, f, exitflag, output, lambda] = om.solve(opt)

%% print results
fprintf('\n-----  METHOD 1 -----');
fprintf('\nf = %g      exitflag = %d\n', f, exitflag);
om.display_soln();

%%-----  METHOD 2  -----
%% assemble model parameters manually
xmin = [ymin; zmin];
xmax = [];
x0 = [y0; z0];

%% solve model
[x, f, exitflag, output, lambda] = qcqps_master(H, c, Q, B, lq, uq, A, b, b, xmin, xmax, x0);
% [x, f, exitflag, output, lambda] = qcqps_master(H, c, Q, B, lq, uq, A, b, b, xmin, xmax, x0, opt)

%% print results
fprintf('\n-----  METHOD 2 -----');
fprintf('\nf = %g      exitflag = %d\n', f, exitflag);
fprintf('\n             var bound shadow prices\n');
fprintf('     x     lambda.lower  lambda.upper\n');
fprintf('%8.4f  %10.4f  %12.4f\n', [x lambda.lower lambda.upper]');
fprintf('\nquadratic constraint shadow prices\n');
fprintf('lambda.mu_lq  lambda.mu_uq\n');
fprintf('%8.4f  %11.4f\n', [lambda.mu_lq lambda.mu_uq]');
fprintf('\nlinear constraint shadow prices\n');
fprintf('lambda.mu_l  lambda.mu_u\n');
fprintf('%8.4f  %11.4f\n', [lambda.mu_l lambda.mu_u]');
