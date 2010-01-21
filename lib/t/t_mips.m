function t_mips(quiet)
%T_MIPS  Tests of MIPS NLP solver

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2010 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

t_begin(51, quiet);

t = 'unconstrained banana function : ';
%% from Matlab Optimization Toolbox's bandem.m
ff = @(x)ff2(x);
function [f, df, d2f] = ff2(x)
    a = 100;
    f = a*(x(2)-x(1)^2)^2+(1-x(1))^2;
    df = [  4*a*(x(1)^3 - x(1)*x(2)) + 2*x(1)-2;
            2*a*(x(2) - x(1)^2)                     ];
    d2f = 4*a*[ 3*x(1)^2 - x(2) + 1/(2*a),  -x(1);
                -x(1)                       1/2       ];
end
x0 = [-1.9; 2];
[x, f, s, out, lam] = mips(ff, x0);
% [x, f, s, out, lam] = mips(ff, x0, [], [], [], [], [], [], [], struct('verbose', 2));
t_ok(s, [t 'success']);
t_is(x, [1; 1], 13, [t 'x']);
t_is(f, 0, 13, [t 'f']);
t_is(out.hist(end).compcond, 0, 6, [t 'compcond']);
t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
t_is(lam.lower, zeros(size(x)), 13, [t 'lam.lower']);
t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

t = 'unconstrained 3-d quadratic : ';
%% from http://www.akiti.ca/QuadProgEx0Constr.html
ff = @(x)ff3(x);
function [f, df, d2f] = ff3(x)
    H = [5 -2 -1; -2 4 3; -1 3 5];
    c = [2; -35; -47];
    f = 1/2 * x'*H*x + c'*x + 5;
    df = H*x + c;
    d2f = H;
end
x0 = [0; 0; 0];
% [x, f, s, out, lam] = mips(ff, x0, [], [], [], [], [], [], [], struct('verbose', 2));
[x, f, s, out, lam] = mips(ff, x0);
t_ok(s, [t 'success']);
t_is(x, [3; 5; 7], 13, [t 'x']);
t_is(f, -244, 13, [t 'f']);
t_is(out.hist(end).compcond, 0, 6, [t 'compcond']);
t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
t_is(lam.lower, zeros(size(x)), 13, [t 'lam.lower']);
t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

t = 'constrained 4-d QP : ';
%% from http://www.uc.edu/sashtml/iml/chap8/sect12.htm
ff = @(x)ff4(x);
function [f, df, d2f] = ff4(x)
    H = [   1003.1  4.3     6.3     5.9;
            4.3     2.2     2.1     3.9;
            6.3     2.1     3.5     4.8;
            5.9     3.9     4.8     10  ];
    c = zeros(4,1);
    f = 1/2 * x'*H*x + c'*x;
    df = H*x + c;
    d2f = H;
end
x0 = [1; 0; 0; 1];
A = [   1       1       1       1;
        0.17    0.11    0.10    0.18    ];
l = [1; 0.10];
u = [1; Inf];
xmin = zeros(4,1);
% [x, f, s, out, lam] = mips(ff, x0, A, l, u, xmin, [], [], [], struct('verbose', 2));
[x, f, s, out, lam] = mips(ff, x0, A, l, u, xmin);
t_ok(s, [t 'success']);
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
% t_ok(s, [t 'success']);
% t_is(x, [0; 2.8; 0.2; 0]/3, 6, [t 'x']);
% t_is(f, 3.29/3, 6, [t 'f']);
% t_is(lam.eqlin, -6.58/3, 6, [t 'lam.eqlin']);
% t_is(lam.ineqlin, 0, 13, [t 'lam.ineqlin']);
% t_is(lam.lower, [2.24;0;0;1.7667], 4, [t 'lam.lower']);
% t_is(lam.upper, [0;0;0;0], 13, [t 'lam.upper']);

t = 'constrained 2-d non-linear : ';
%% from http://en.wikipedia.org/wiki/Nonlinear_programming#2-dimensional_example
ff = @(x)ff5(x);
gh = @(x)gh5(x);
hess = @(x, lam)hess5(x, lam);
function [f, df, d2f] = ff5(x)
    c = -[1; 1];
    f = c'*x;
    df = c;
    d2f = zeros(2,2);
end
function [g, h, dg, dh] = gh5(x)
    g = [ -1 -1; 1 1] * x.^2 + [1; -2];
    dg = 2 * [-x(1) x(1); -x(2) x(2)];
    h = []; dh = [];
end
function Lxx = hess5(x, lam)
    mu = lam.ineqnonlin;
    Lxx = 2*[-1 1]*mu*eye(2);
end
x0 = [1.1; 0];
xmin = -10 * ones(2, 1);
% xmax = 3 * ones(2, 1);
% [x, f, s, out, lam] = mips(ff, x0, [], [], [], xmin, [], gh, hess, struct('verbose', 2));
[x, f, s, out, lam] = mips(ff, x0, [], [], [], xmin, [], gh, hess);
t_ok(s, [t 'success']);
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
% % [x, f, s, out, lam] = fmincon(ff, x0, [], [], [], [], xmin, [], gh, fmoptions);
% [x, f, s, out, lam] = fmincon(ff, x0, [], [], [], [], [], [], gh);
% t_ok(s, [t 'success']);
% t_is(x, [1; 1], 4, [t 'x']);
% t_is(f, -2, 6, [t 'f']);
% t_is(lam.ineqnonlin, [0;0.5], 6, [t 'lam.ineqnonlin']);

t = 'constrained 3-d non-linear : ';
%% from http://en.wikipedia.org/wiki/Nonlinear_programming#3-dimensional_example
ff = @(x)ff6(x);
gh = @(x)gh6(x);
hess = @(x, lam)hess6(x, lam);
function [f, df, d2f] = ff6(x)
    f = -x(1)*x(2) - x(2)*x(3);
    df = -[x(2); x(1)+x(3); x(2)];
    d2f = -[0 1 0; 1 0 1; 0 1 0];
end
function [g, h, dg, dh] = gh6(x)
    g = [ 1 -1 1; 1 1 1] * x.^2 + [-2; -10];
    dg = 2 * [x(1) x(1); -x(2) x(2); x(3) x(3)];
    h = []; dh = [];
end
function Lxx = hess6(x, lam)
    mu = lam.ineqnonlin;
    Lxx = [2*[1 1]*mu -1 0; 0 2*[-1 1]*mu -1; 0 -1 2*[1 1]*mu];
end
x0 = [1; 1; 0];
% [x, f, s, out, lam] = mips(ff, x0, [], [], [], [], [], gh, hess, struct('verbose', 2, 'comptol', 1e-9));
[x, f, s, out, lam] = mips(ff, x0, [], [], [], [], [], gh, hess);
t_ok(s, [t 'success']);
t_is(x, [1.58113883; 2.23606798; 1.58113883], 6, [t 'x']);
t_is(f, -5*sqrt(2), 6, [t 'f']);
t_is(out.hist(end).compcond, 0, 6, [t 'compcond']);
t_is(lam.ineqnonlin, [0;sqrt(2)/2], 7, [t 'lam.ineqnonlin']);
t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
t_is(lam.lower, zeros(size(x)), 13, [t 'lam.lower']);
t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);
% %% check with fmincon (for dev testing only)
% % fmoptions = optimset('Algorithm', 'interior-point');
% % [x, f, s, out, lam] = fmincon(ff, x0, [], [], [], [], xmin, [], gh, fmoptions);
% [x, f, s, out, lam] = fmincon(ff, x0, [], [], [], [], [], [], gh);
% t_ok(s, [t 'success']);
% t_is(x, [1.58113883; 2.23606798; 1.58113883], 4, [t 'x']);
% t_is(f, -5*sqrt(2), 8, [t 'f']);
% t_is(lam.ineqnonlin, [0;sqrt(2)/2], 8, [t 'lam.ineqnonlin']);

t = 'constrained 3-d non-linear (struct) : ';
p = struct('f', ff, 'x0', x0, 'gh', gh, 'hess', hess);
[x, f, s, out, lam] = mips(p);
t_ok(s, [t 'success']);
t_is(x, [1.58113883; 2.23606798; 1.58113883], 6, [t 'x']);
t_is(f, -5*sqrt(2), 6, [t 'f']);
t_is(out.hist(end).compcond, 0, 6, [t 'compcond']);
t_is(lam.ineqnonlin, [0;sqrt(2)/2], 7, [t 'lam.ineqnonlin']);
t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
t_is(lam.lower, zeros(size(x)), 13, [t 'lam.lower']);
t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

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
% [x, f, s, out, lam] = fmincon(ff, x0, [-A; A], [-l; u], [], [], [], [], [], fmoptions);
% t_is(x, [24; 12; 12], 13, t);
% t_is(f, -3456, 13, t);


end
