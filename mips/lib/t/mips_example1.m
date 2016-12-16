f_fcn = @(x)banana(x, 100);
x0 = [-1.9; 2];
[x, f] = mips(f_fcn, x0)
