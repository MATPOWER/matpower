%% variable initial values
y0 = [1; 0];
z0 = [0; 1];

%% variable lower bounds
ymin = [0; 0];
zmin = [0; 0];

%% constraint data
A1 = [ 1 1 1 1 ];               b1 = 1;
A2 = [ 0.17 0.11 0.10 0.18 ];   l2 = 0.1;

%% quadratic cost coefficients
Q = [   1003.1 4.3 6.3 5.9;
        4.3 2.2 2.1 3.9;
        6.3 2.1 3.5 4.8;
        5.9 3.9 4.8 10  ];

%% solver options
opt = struct('verbose', 2, 'alg', 'MIPS');

%%-----  METHOD 1  -----
%% build model
om = opt_model;
om.add_var('y', 2, y0, ymin);
om.add_var('z', 2, z0, zmin);
om.add_lin_constraint('lincon1', A1, b1, b1, {'y', 'z'});
om.add_lin_constraint('lincon2', A2, l2, [], {'y', 'z'});
om.add_quad_cost('cost', Q, [], [], {'y', 'z'});

%% solve model
[x, f, exitflag, output, lambda] = om.solve();
% [x, f, exitflag, output, lambda] = om.solve(opt)

%% print results
fprintf('\n-----  METHOD 1 -----');
fprintf('\nf = %g   exitflag = %d\n', f, exitflag);
fprintf('\nx = \n');
fprintf('   %.4f\n', x);
fprintf('\nlambda.lower (var bound shadow price) =\n');
fprintf('   %.4f\n', lambda.lower);
fprintf('\nlambda.mu_l (constraint shadow price) =\n');
fprintf('   %.4f\n', lambda.mu_l);

%%-----  METHOD 2  -----
%% assemble model parameters manually
xmin = [ymin; zmin];
x0 = [y0; z0];
A = [ A1; A2 ];
l = [ b1; l2 ];
u = [ b1; Inf ];

%% solve model
[x, f, exitflag, output, lambda] = qps_master(Q, [], A, l, u, xmin, [], x0);
% [x, f, exitflag, output, lambda] = qps_master(Q, [], A, l, u, xmin, [], x0, opt)

%% print results
fprintf('\n-----  METHOD 2 -----');
fprintf('\nf = %g   exitflag = %d\n', f, exitflag);
fprintf('\nx = \n');
fprintf('   %.4f\n', x);
fprintf('\nlambda.lower (var bound shadow price) =\n');
fprintf('   %.4f\n', lambda.lower);
fprintf('\nlambda.mu_l (constraint shadow price) =\n');
fprintf('   %.4f\n', lambda.mu_l);
