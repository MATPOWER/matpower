function t_mips_pardiso(quiet)
%T_MIPS_PARDISO  Tests of MIPS NLP solver, using PARDISO as linear solver.

%   MIPS
%   Copyright (c) 2010-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MIPS.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mips for more info.

if nargin < 1
    quiet = 0;
end

num_tests = 60;

t_begin(num_tests, quiet);

if have_pardiso()

opt = struct('linsolver', 'PARDISO');

t = 'unconstrained banana function : ';
%% from MATLAB Optimization Toolbox's bandem.m
f_fcn = @(x)f2(x);
x0 = [-1.9; 2];
% [x, f, s, out, lam] = mips(f_fcn, x0, [], [], [], [], [], [], [], struct('verbose', 2));
[x, f, s, out, lam] = mips(f_fcn, x0, [], [], [], [], [], [], [], opt);
t_is(s, 1, 13, [t 'success']);
t_is(x, [1; 1], 13, [t 'x']);
t_is(f, 0, 13, [t 'f']);
t_is(out.hist(end).compcond, 0, 6, [t 'compcond']);
t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
t_is(lam.lower, zeros(size(x)), 13, [t 'lam.lower']);
t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

t = 'unconstrained 3-d quadratic : ';
%% from http://www.akiti.ca/QuadProgEx0Constr.html
f_fcn = @(x)f3(x);
x0 = [0; 0; 0];
% [x, f, s, out, lam] = mips(f_fcn, x0, [], [], [], [], [], [], [], struct('verbose', 2));
[x, f, s, out, lam] = mips(f_fcn, x0, [], [], [], [], [], [], [], opt);
t_is(s, 1, 13, [t 'success']);
t_is(x, [3; 5; 7], 13, [t 'x']);
t_is(f, -244, 13, [t 'f']);
t_is(out.hist(end).compcond, 0, 6, [t 'compcond']);
t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
t_is(lam.lower, zeros(size(x)), 13, [t 'lam.lower']);
t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

t = 'constrained 4-d QP : ';
%% from http://www.jmu.edu/docs/sasdoc/sashtml/iml/chap8/sect12.htm
f_fcn = @(x)f4(x);
x0 = [1; 0; 0; 1];
A = [   1       1       1       1;
        0.17    0.11    0.10    0.18    ];
l = [1; 0.10];
u = [1; Inf];
xmin = zeros(4,1);
% [x, f, s, out, lam] = mips(f_fcn, x0, A, l, u, xmin, [], [], [], struct('verbose', 2));
[x, f, s, out, lam] = mips(f_fcn, x0, A, l, u, xmin, [], [], [], opt);
t_is(s, 1, 13, [t 'success']);
t_is(x, [0; 2.8; 0.2; 0]/3, 6, [t 'x']);
t_is(f, 3.29/3, 6, [t 'f']);
t_is(out.hist(end).compcond, 0, 6, [t 'compcond']);
t_is(lam.mu_l, [6.58;0]/3, 6, [t 'lam.mu_l']);
t_is(lam.mu_u, [0;0], 13, [t 'lam.mu_u']);
t_is(lam.lower, [2.24;0;0;1.7667], 4, [t 'lam.lower']);
t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

H = [   1003.1  4.3     6.3     5.9;
        4.3     2.2     2.1     3.9;
        6.3     2.1     3.5     4.8;
        5.9     3.9     4.8     10  ];
c = zeros(4,1);
% %% check with quadprog (for dev testing only)
% [x, f, s, out, lam] = quadprog(H,c,-A(2,:), -0.10, A(1,:), 1, xmin);
% t_is(s, 1, 13, [t 'success']);
% t_is(x, [0; 2.8; 0.2; 0]/3, 6, [t 'x']);
% t_is(f, 3.29/3, 6, [t 'f']);
% t_is(lam.eqlin, -6.58/3, 6, [t 'lam.eqlin']);
% t_is(lam.ineqlin, 0, 13, [t 'lam.ineqlin']);
% t_is(lam.lower, [2.24;0;0;1.7667], 4, [t 'lam.lower']);
% t_is(lam.upper, [0;0;0;0], 13, [t 'lam.upper']);

t = 'constrained 2-d nonlinear : ';
%% from http://en.wikipedia.org/wiki/Nonlinear_programming#2-dimensional_example
f_fcn = @(x)f5(x);
gh_fcn = @(x)gh5(x);
hess_fcn = @(x, lam, cost_mult)hess5(x, lam, cost_mult);
x0 = [1.1; 0];
xmin = zeros(2, 1);
% xmax = 3 * ones(2, 1);
% [x, f, s, out, lam] = mips(f_fcn, x0, [], [], [], xmin, [], gh_fcn, hess_fcn, struct('verbose', 2));
[x, f, s, out, lam] = mips(f_fcn, x0, [], [], [], xmin, [], gh_fcn, hess_fcn, opt);
t_is(s, 1, 13, [t 'success']);
t_is(x, [1; 1], 6, [t 'x']);
t_is(f, -2, 6, [t 'f']);
t_is(out.hist(end).compcond, 0, 6, [t 'compcond']);
t_is(lam.ineqnonlin, [0;0.5], 6, [t 'lam.ineqnonlin']);
t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
t_is(lam.lower, zeros(size(x)), 13, [t 'lam.lower']);
t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);
% %% check with fmincon (for dev testing only)
% % fmoptions = optimset('Algorithm', 'interior-point');
% % [x, f, s, out, lam] = fmincon(f_fcn, x0, [], [], [], [], xmin, [], gh_fcn, fmoptions);
% [x, f, s, out, lam] = fmincon(f_fcn, x0, [], [], [], [], [], [], gh_fcn);
% t_is(s, 1, 13, [t 'success']);
% t_is(x, [1; 1], 4, [t 'x']);
% t_is(f, -2, 6, [t 'f']);
% t_is(lam.ineqnonlin, [0;0.5], 6, [t 'lam.ineqnonlin']);

t = 'constrained 3-d nonlinear : ';
%% from http://en.wikipedia.org/wiki/Nonlinear_programming#3-dimensional_example
f_fcn = @(x)f6(x);
gh_fcn = @(x)gh6(x);
hess_fcn = @(x, lam, cost_mult)hess6(x, lam, cost_mult);
x0 = [1; 1; 0];
% [x, f, s, out, lam] = mips(f_fcn, x0, [], [], [], [], [], gh_fcn, hess_fcn, struct('verbose', 2, 'comptol', 1e-9));
[x, f, s, out, lam] = mips(f_fcn, x0, [], [], [], [], [], gh_fcn, hess_fcn, opt);
t_is(s, 1, 13, [t 'success']);
t_is(x, [1.58113883; 2.23606798; 1.58113883], 6, [t 'x']);
t_is(f, -5*sqrt(2), 6, [t 'f']);
t_is(out.hist(end).compcond, 0, 6, [t 'compcond']);
t_is(lam.ineqnonlin, [0;sqrt(2)/2], 7, [t 'lam.ineqnonlin']);
t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
t_is(lam.lower, zeros(size(x)), 13, [t 'lam.lower']);
t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

t = 'constrained 3-d nonlinear (struct) : ';
p = struct('f_fcn', f_fcn, 'x0', x0, 'gh_fcn', gh_fcn, 'hess_fcn', hess_fcn);
[x, f, s, out, lam] = mips(p);
t_is(s, 1, 13, [t 'success']);
t_is(x, [1.58113883; 2.23606798; 1.58113883], 6, [t 'x']);
t_is(f, -5*sqrt(2), 6, [t 'f']);
t_is(out.hist(end).compcond, 0, 6, [t 'compcond']);
t_is(lam.ineqnonlin, [0;sqrt(2)/2], 7, [t 'lam.ineqnonlin']);
t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
t_is(lam.lower, zeros(size(x)), 13, [t 'lam.lower']);
t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

t = 'constrained 4-d nonlinear : ';
%% Hock & Schittkowski test problem #71
f_fcn = @(x)f7(x);
gh_fcn = @(x)gh7(x);
hess_fcn = @(x, lam, sigma)hess7(x, lam, sigma);
x0 = [1; 5; 5; 1];
xmin = ones(4, 1);
xmax = 5 * xmin;
% [x, f, s, out, lam] = mips(f_fcn, x0, [], [], [], xmin, xmax, gh_fcn, hess_fcn, struct('verbose', 2, 'comptol', 1e-9));
[x, f, s, out, lam] = mips(f_fcn, x0, [], [], [], xmin, xmax, gh_fcn, hess_fcn, opt);
t_is(s, 1, 13, [t 'success']);
t_is(x, [1; 4.7429994; 3.8211503; 1.3794082], 6, [t 'x']);
t_is(f, 17.0140173, 6, [t 'f']);
t_is(lam.eqnonlin, 0.1614686, 5, [t 'lam.eqnonlin']);
t_is(lam.ineqnonlin, 0.55229366, 5, [t 'lam.ineqnonlin']);
t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
t_is(lam.lower, [1.08787121024; 0; 0; 0], 5, [t 'lam.lower']);
t_is(lam.upper, zeros(size(x)), 7, [t 'lam.upper']);

else
    t_skip(num_tests, 'PARDISO not available.');
end

t_end;


% %%-----  eg99 : linearly constrained fmincon example, mips can't solve  -----
% function [f, df, d2f] = eg99(x)
% f = -x(1)*x(2)*x(3);
% df = -[ x(2)*x(3);
%         x(1)*x(3);
%         x(1)*x(2)   ];
% d2f = -[    0       x(3)    x(2);
%             x(3)    0       x(1);
%             x(2)    x(1)    0   ];
% end
% 
% x0 = [10;10;10];
% A = [1 2 2];
% l = 0;
% u = 72;
% fmoptions = optimset('Display', 'testing');
% fmoptions = optimset(fmoptions, 'Algorithm', 'interior-point');
% [x, f, s, out, lam] = fmincon(f_fcn, x0, [-A; A], [-l; u], [], [], [], [], [], fmoptions);
% t_is(x, [24; 12; 12], 13, t);
% t_is(f, -3456, 13, t);


%% unconstrained banana function
%% from MATLAB Optimization Toolbox's bandem.m
function [f, df, d2f] = f2(x)
    a = 100;
    f = a*(x(2)-x(1)^2)^2+(1-x(1))^2;
    df = [  4*a*(x(1)^3 - x(1)*x(2)) + 2*x(1)-2;
            2*a*(x(2) - x(1)^2)                     ];
    d2f = 4*a*[ 3*x(1)^2 - x(2) + 1/(2*a),  -x(1);
                -x(1)                       1/2       ];


%% unconstrained 3-d quadratic
%% from http://www.akiti.ca/QuadProgEx0Constr.html
function [f, df, d2f] = f3(x)
    H = [5 -2 -1; -2 4 3; -1 3 5];
    c = [2; -35; -47];
    f = 1/2 * x'*H*x + c'*x + 5;
    df = H*x + c;
    d2f = H;


%% constrained 4-d QP
%% from http://www.jmu.edu/docs/sasdoc/sashtml/iml/chap8/sect12.htm
function [f, df, d2f] = f4(x)
    H = [   1003.1  4.3     6.3     5.9;
            4.3     2.2     2.1     3.9;
            6.3     2.1     3.5     4.8;
            5.9     3.9     4.8     10  ];
    c = zeros(4,1);
    f = 1/2 * x'*H*x + c'*x;
    df = H*x + c;
    d2f = H;


%% constrained 2-d nonlinear
%% from http://en.wikipedia.org/wiki/Nonlinear_programming#2-dimensional_example
function [f, df, d2f] = f5(x)
    c = -[1; 1];
    f = c'*x;
    df = c;
    d2f = zeros(2,2);

function [h, g, dh, dg] = gh5(x)
    h = [ -1 -1; 1 1] * x.^2 + [1; -2];
    dh = 2 * [-x(1) x(1); -x(2) x(2)];
    g = []; dg = [];

function Lxx = hess5(x, lam, cost_mult)
    mu = lam.ineqnonlin;
    Lxx = 2*[-1 1]*mu*eye(2);


%% constrained 3-d nonlinear
%% from http://en.wikipedia.org/wiki/Nonlinear_programming#3-dimensional_example
function [f, df, d2f] = f6(x)
    f = -x(1)*x(2) - x(2)*x(3);
    df = -[x(2); x(1)+x(3); x(2)];
    d2f = -[0 1 0; 1 0 1; 0 1 0];

function [h, g, dh, dg] = gh6(x)
    h = [ 1 -1 1; 1 1 1] * x.^2 + [-2; -10];
    dh = 2 * [x(1) x(1); -x(2) x(2); x(3) x(3)];
    g = []; dg = [];

function Lxx = hess6(x, lam, cost_mult)
    if nargin < 3, cost_mult = 1; end
    mu = lam.ineqnonlin;
    Lxx = cost_mult * [0 -1 0; -1 0 -1; 0 -1 0] + ...
            [2*[1 1]*mu 0 0; 0 2*[-1 1]*mu 0; 0 0 2*[1 1]*mu];


%% constrained 4-d nonlinear
%% Hock & Schittkowski test problem #71
function [f, df, d2f] = f7(x)
    f = x(1)*x(4)*sum(x(1:3)) + x(3);
    df = [ x(1)*x(4) + x(4)*sum(x(1:3));
           x(1)*x(4);
           x(1)*x(4) + 1;
           x(1)*sum(x(1:3)) ];
    d2f = sparse([ 2*x(4)        x(4)   x(4)  2*x(1)+x(2)+x(3);
              x(4)               0      0     x(1);
              x(4)               0      0     x(1);
              2*x(1)+x(2)+x(3)  x(1)  x(1)    0
        ]);

function [h, g, dh, dg] = gh7(x)
    g = sum(x.^2) - 40;
    h = -prod(x) + 25;
    dg = 2*x;
    dh = -prod(x)./x;

function Lxx = hess7(x, lam, sigma)
    if nargin < 3, sigma = 1; end
    lambda = lam.eqnonlin;
    mu     = lam.ineqnonlin;
    [f, df, d2f] = f7(x);
    Lxx = sigma * d2f + lambda*2*speye(4) - ...
       mu*sparse([      0     x(3)*x(4) x(2)*x(4) x(2)*x(3);
                    x(3)*x(4)     0     x(1)*x(4) x(1)*x(3);
                    x(2)*x(4) x(1)*x(4)     0     x(1)*x(2);
                    x(2)*x(3) x(1)*x(3) x(1)*x(2)     0  ]);


function TorF = have_pardiso()
TorF = have_pardiso_object() || have_pardiso_legacy();

function TorF = have_pardiso_object()
TorF = exist('pardiso', 'file') == 2;
if TorF
    try
        id = 1;
        A = sparse([1 2; 3 4]);
        b = [1;1];
        p = pardiso(id, 1, 0);
        p.factorize(id, A);
        x = p.solve(id, A, b);
        p.free(id);
        p.clear();
        if any(x ~= [-1; 1])
            TorF = 0;
        end
    catch
        TorF = 0;
    end
end

function TorF = have_pardiso_legacy()
TorF = exist('pardisoinit', 'file') == 3 && ...
        exist('pardisoreorder', 'file') == 3 && ...
        exist('pardisofactor', 'file') == 3 && ...
        exist('pardisosolve', 'file') == 3 && ...
        exist('pardisofree', 'file') == 3;
if TorF
    try
        A = sparse([1 2; 3 4]);
        b = [1;1];
        info = pardisoinit(11, 0);
        info = pardisoreorder(A, info, false);
        info = pardisofactor(A, info, false);
        [x, info] = pardisosolve(A, b, info, false);
        pardisofree(info);
        if any(x ~= [-1; 1])
            TorF = 0;
        end
    catch
        info
        TorF = 0;
    end
end
