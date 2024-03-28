function t_nlps_master(quiet)
% t_nlps_master - Tests of NLP solvers via nlps_master.

%   MP-Opt-Model
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin < 1
    quiet = 0;
end

%%  alg         name        check           opts
cfg = {
    {'DEFAULT', 'DEFAULT',  [],             []                          },
    {'MIPS',    'MIPS',     [],             []                          },
    {'MIPS',    'sc-MIPS',  [],             struct('step_control', 1)   },
    {'FMINCON', 'fmincon-2','fmincon_ipm',  struct('alg', 2)            },
    {'FMINCON', 'fmincon-3','fmincon_ipm',  struct('alg', 3)            },
    {'FMINCON', 'fmincon-4','fmincon_ipm',  struct('alg', 4)            },
    {'FMINCON', 'fmincon-5','fmincon_ipm',  struct('alg', 5)            },
    {'FMINCON', 'fmincon-6','fmincon_ipm',  struct('alg', 6)            },
    {'IPOPT',   'IPOPT',    'ipopt',        []                          },
    {'KNITRO',  'Knitro',   'knitro',       []                          },
};
if have_feature('matlab') && have_feature('matlab', 'vnum') <= 8.003
    cfg([5;8]) = [];    %% Mac MATLAB 7.14-8.2 crash w/ fmincon alg 3, fails w/6
end
% cfg = {
%     {'KNITRO',  'Knitro',   'knitro',       []                          },
% };

n = 54;

t_begin(n*length(cfg), quiet);

for k = 1:length(cfg)
    alg   = cfg{k}{1};
    name  = cfg{k}{2};
    check = cfg{k}{3};
    opts  = cfg{k}{4};
    if ~isempty(check) && ~have_feature(check)
        t_skip(n, sprintf('%s not installed', name));
    else
        opt = struct('verbose', 0, 'alg', alg);
        switch alg
            case {'DEFAULT', 'MIPS'}
                opt.mips_opt = opts;
                opt.mips_opt.comptol = 1e-8;
            case 'FMINCON'
                opt.fmincon_opt = opts;
%                 opt.fmincon_opt.opts.OptimalityTolerance = 1e-10;
%                 opt.fmincon_opt.opts.TolCon = 1e-8;
%                 opt.fmincon_opt.tol_x = 1e-8;
%                 opt.fmincon_opt.opts.MaxFunctionEvaluations = 10000;
%                 opt.fmincon_opt.opts.FiniteDifferenceType = 'central';
                opt.fmincon_opt.tol_f = 1e-9;
            case 'IPOPT'
                opt.ipopt_opt = opts;
%                 opt.ipopt_opt.acceptable_tol = 1e-11;
%                 opt.ipopt_opt.acceptable_compl_inf_tol = 1e-11;
%                 opt.ipopt_opt.acceptable_constr_viol_tol = 1e-11;
%                 opt.ipopt_opt.dual_inf_tol = 1e-9;
                opt.ipopt_opt.compl_inf_tol   = 1e-9;
%                 opt.ipopt_opt.acceptable_obj_change_tol   = 1e-10;
%                 opt.ipopt_opt.linear_solver   = 'ma57';
            case 'KNITRO'
                opt.knitro_opt = opts;
                opt.knitro_opt.tol_f = 1e-8;
        end

t = sprintf('%s - unconstrained banana function : ', name);
%% from MATLAB Optimization Toolbox's bandem.m
f_fcn = @(x)f2(x);
x0 = [-1.9; 2];
[x, f, s, out, lam] = nlps_master(f_fcn, x0, [], [], [], [], [], [], [], opt);
t_ok(s > 0, [t 'success']);
t_is(x, [1; 1], 8, [t 'x']);
t_is(f, 0, 13, [t 'f']);
t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
t_is(lam.lower, zeros(size(x)), 13, [t 'lam.lower']);
t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

t = sprintf('%s - unconstrained 3-d quadratic : ', name);
%% from http://www.akiti.ca/QuadProgEx0Constr.html
f_fcn = @(x)f3(x);
x0 = [0; 0; 0];
[x, f, s, out, lam] = nlps_master(f_fcn, x0, [], [], [], [], [], [], [], opt);
t_ok(s > 0, [t 'success']);
t_is(x, [3; 5; 7], 6, [t 'x']);
t_is(f, -244, 12, [t 'f']);
t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
t_is(lam.lower, zeros(size(x)), 13, [t 'lam.lower']);
t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

t = sprintf('%s - constrained 4-d QP : ', name);
%% from https://v8doc.sas.com/sashtml/iml/chap8/sect12.htm
f_fcn = @(x)f4(x);
x0 = [1; 0; 0; 1];
A = [   1       1       1       1;
        0.17    0.11    0.10    0.18    ];
l = [1; 0.10];
u = [1; Inf];
xmin = zeros(4,1);
[x, f, s, out, lam] = nlps_master(f_fcn, x0, A, l, u, xmin, [], [], [], opt);
t_ok(s > 0, [t 'success']);
t_is(x, [0; 2.8; 0.2; 0]/3, 6, [t 'x']);
t_is(f, 3.29/3, 6, [t 'f']);
t_is(lam.mu_l, [6.58;0]/3, 5, [t 'lam.mu_l']);
t_is(lam.mu_u, [0;0], 13, [t 'lam.mu_u']);
t_is(lam.lower, [2.24;0;0;1.7667], 4, [t 'lam.lower']);
t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

t = sprintf('%s - constrained 2-d nonlinear : ', name);
%% from https://en.wikipedia.org/wiki/Nonlinear_programming#2-dimensional_example
f_fcn = @(x)f5(x);
gh_fcn = @(x)gh5(x);
hess_fcn = @(x, lam, cost_mult)hess5(x, lam, cost_mult);
x0 = [1.1; 0];
xmin = zeros(2, 1);
% xmax = 3 * ones(2, 1);
[x, f, s, out, lam] = nlps_master(f_fcn, x0, [], [], [], xmin, [], gh_fcn, hess_fcn, opt);
t_ok(s > 0, [t 'success']);
t_is(x, [1; 1], 6, [t 'x']);
t_is(f, -2, 6, [t 'f']);
t_is(lam.ineqnonlin, [0;0.5], 6, [t 'lam.ineqnonlin']);
t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
t_is(lam.lower, zeros(size(x)), 9, [t 'lam.lower']);
t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

t = sprintf('%s - constrained 3-d nonlinear : ', name);
if strcmp(name, 'fmincon-5') && have_feature('matlab', 'vnum') >= 9.011
    t_skip(16, [t 'known issue w/R2021b']);
else
    %% from https://en.wikipedia.org/wiki/Nonlinear_programming#3-dimensional_example
    f_fcn = @(x)f6(x);
    gh_fcn = @(x)gh6(x);
    hess_fcn = @(x, lam, cost_mult)hess6(x, lam, cost_mult);
    x0 = [1; 1; 0];
    [x, f, s, out, lam] = nlps_master(f_fcn, x0, [], [], [], [], [], gh_fcn, hess_fcn, opt);
    t_ok(s > 0, [t 'success']);
    t_is(x, [1.58113883; 2.23606798; 1.58113883], 6, [t 'x']);
    t_is(f, -5*sqrt(2), 6, [t 'f']);
    t_is(lam.ineqnonlin, [0;sqrt(2)/2], 7, [t 'lam.ineqnonlin']);
    t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
    t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
    t_is(lam.lower, zeros(size(x)), 13, [t 'lam.lower']);
    t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

    t = sprintf('%s - constrained 3-d nonlinear (struct) : ', name);
    p = struct('f_fcn', f_fcn, 'x0', x0, 'gh_fcn', gh_fcn, 'hess_fcn', hess_fcn, 'opt', opt);
    [x, f, s, out, lam] = nlps_master(p);
    t_ok(s > 0, [t 'success']);
    t_is(x, [1.58113883; 2.23606798; 1.58113883], 6, [t 'x']);
    t_is(f, -5*sqrt(2), 6, [t 'f']);
    t_is(lam.ineqnonlin, [0;sqrt(2)/2], 7, [t 'lam.ineqnonlin']);
    t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
    t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
    t_is(lam.lower, zeros(size(x)), 13, [t 'lam.lower']);
    t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);
end

t = sprintf('%s - constrained 4-d nonlinear : ', name);
%% Hock & Schittkowski test problem #71
f_fcn = @(x)f7(x);
gh_fcn = @(x)gh7(x);
hess_fcn = @(x, lam, sigma)hess7(x, lam, sigma);
x0 = [1; 5; 5; 1];
xmin = ones(4, 1);
xmax = 5 * xmin;
[x, f, s, out, lam] = nlps_master(f_fcn, x0, [], [], [], xmin, xmax, gh_fcn, hess_fcn, opt);
t_ok(s > 0, [t 'success']);
t_is(x, [1; 4.7429994; 3.8211503; 1.3794082], 6, [t 'x']);
t_is(f, 17.0140173, 6, [t 'f']);
t_is(lam.eqnonlin, 0.1614686, 5, [t 'lam.eqnonlin']);
t_is(lam.ineqnonlin, 0.55229366, 5, [t 'lam.ineqnonlin']);
t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
t_is(lam.lower, [1.08787121024; 0; 0; 0], 5, [t 'lam.lower']);
t_is(lam.upper, zeros(size(x)), 7, [t 'lam.upper']);
    end
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
%% from https://v8doc.sas.com/sashtml/iml/chap8/sect12.htm
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
%% from https://en.wikipedia.org/wiki/Nonlinear_programming#2-dimensional_example
function [f, df, d2f] = f5(x)
    c = -[1; 1];
    f = c'*x;
    df = c;
    d2f = zeros(2,2);

function [h, g, dh, dg] = gh5(x)
    [h, dh] = h5(x);
    dh = dh';
    g = []; dg = [];

function [h, dh] = h5(x)
    h = [ -1 -1; 1 1] * x.^2 + [1; -2];
    dh = 2 * [-x(1) -x(2); x(1) x(2)];

function Lxx = hess5(x, lam, cost_mult)
    Lxx = hess5a(x, lam.ineqnonlin);

function Lxx = hess5a(x, mu)
    Lxx = 2*[-1 1]*mu*eye(2);

%% constrained 3-d nonlinear
%% from https://en.wikipedia.org/wiki/Nonlinear_programming#3-dimensional_example
function [f, df, d2f] = f6(x)
    f = -x(1)*x(2) - x(2)*x(3);
    df = -[x(2); x(1)+x(3); x(2)];
    d2f = -[0 1 0; 1 0 1; 0 1 0];

function [h, g, dh, dg] = gh6(x)
    [h, dh] = h6(x);
    dh = dh';
    g = []; dg = [];

function [h, dh] = h6(x)
    h = [ 1 -1 1; 1 1 1] * x.^2 + [-2; -10];
    dh = 2 * [x(1) -x(2) x(3); x(1) x(2) x(3)];

function Lxx = hess6(x, lam, cost_mult)
    if nargin < 3, cost_mult = 1; end
    Lxx = cost_mult * [0 -1 0; -1 0 -1; 0 -1 0] + hess6a(x, lam.ineqnonlin);

function d2h = hess6a(x, mu)
    d2h = [2*[1 1]*mu 0 0; 0 2*[-1 1]*mu 0; 0 0 2*[1 1]*mu];


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
    [g, dg] = g7(x);
    [h, dh] = h7(x);
    dg = dg';
    dh = dh';

function [g, dg] = g7(x)
    g = sum(x.^2) - 40;
    dg = 2*x';

function [h, dh] = h7(x)
    h = -prod(x) + 25;
    dh = -(prod(x)./x)';

function Lxx = hess7(x, lam, sigma)
    if nargin < 3, sigma = 1; end
    [f, df, d2f] = f7(x);
    Lxx = sigma * d2f + hess7g(x, lam.eqnonlin) + hess7h(x, lam.ineqnonlin);

function d2g = hess7g(x, lam)
    d2g = lam*2*speye(4);

function d2h = hess7h(x, mu)
    d2h = -mu*sparse([      0     x(3)*x(4) x(2)*x(4) x(2)*x(3);
                        x(3)*x(4)     0     x(1)*x(4) x(1)*x(3);
                        x(2)*x(4) x(1)*x(4)     0     x(1)*x(2);
                        x(2)*x(3) x(1)*x(3) x(1)*x(2)     0  ]);
