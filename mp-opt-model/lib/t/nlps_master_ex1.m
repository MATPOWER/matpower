function nlps_master_ex1
% nlps_master_ex1 - Example of unconstrained nonlinear optimization.
%
% Minimizes MATLAB's :func:`banana` function using nlps_master.
f_fcn = @(x)banana(x, 100);
x0 = [-1.9; 2];
[x, f] = nlps_master(f_fcn, x0)
