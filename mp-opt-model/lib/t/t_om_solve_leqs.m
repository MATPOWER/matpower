function t_om_solve_leqs(quiet)
%T_OM_SOLVE_LEQS  Tests of LEQ solvers via OM.SOLVE().

%   MP-Opt-Model
%   Copyright (c) 2010-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin < 1
    quiet = 0;
end

%%  alg     name        check   opts
cfg = {
    {'',    'default',  []      []  },
    {'\',   '\',        []      []  },
    {'LU',  'LU',       []      []  },
%     {'PARDISO',  'PARDISO',       []      []  },
};

A1 = [2 1; -3 5];
b1 = [8; 1];
x1 = [3;2];

A2 = sparse([2 -1 0; -3 1 -2; 0 5 -4]);
b2 = [-5; 1; -7];
x2 = [-2; 1; 3];

n = 7;

t_begin(10+n*length(cfg), quiet);

for k = 1:length(cfg)
    alg   = cfg{k}{1};
    name  = cfg{k}{2};

    t = sprintf('%s - full A matrix : ', name);
    om = opt_model;
    om.add_var('x', 2, zeros(size(x1)));
    om.add_lin_constraint('A', A1, b1, b1);
    if isempty(alg)
        opt = struct();
    else
        opt = struct('leq_opt', struct('solver', alg));
    end
    x = om.solve(opt);
    t_is(x, x1, 14, [t 'x']);
    t_ok(~isfield(om.soln, 'var'), [t 'no parse_soln() outputs']);
    
    t = sprintf('%s - sparse A matrix : ', name);
    om = opt_model;
    om.add_var('x', 3, zeros(size(x2)));
    om.add_lin_constraint('A12', A2(1:2, :), b2(1:2), b2(1:2));
    om.add_lin_constraint('A3', A2(3, :), b2(3), b2(3));
    opt = struct('leq_opt', struct('solver', alg), 'parse_soln', 1);
    [x, f, e, out, jac] = om.solve(opt);
    t_is(x, x2, 14, [t 'x']);
    t_is(f, A2*x-b2, 14, [t 'f']);
    t_is(e, 1, 14, [t 'exitflag']);
    t_ok(strcmp(out.alg, alg), [t 'output']);
    t_is(jac, A2, 14, [t 'jac']);
    opt.parse_soln = 0;
end

t = 'om.soln.';
t_is(om.soln.x, x, 14, [t 'x']);
t_is(om.soln.f, f, 14, [t 'f']);
t_is(om.soln.eflag, e, 14, [t 'eflag']);
t_ok(strcmp(om.soln.output.alg, out.alg), [t 'output.alg']);
t_is(om.soln.jac, jac, 14, [t 'jac']);

t = 'om.get_soln(''var'', ''x'') : ';
t_is(om.get_soln('var', 'x'), x, 14, [t 'x']);

t = 'om.get_soln(''lin'', ''A12'') : ';
f12 = om.get_soln('lin', 'A12');
t_is(f12, f(1:2), 14, [t 'A12 * x - u12']);

t = 'om.get_soln(''lin'', ''f'', ''A3'') : ';
g = om.get_soln('lin', 'Ax_u', 'A3');
t_is(g, f(3), 14, [t 'A3 * x - u3']);

t = 'om.get_soln(''lin'', ''l_Ax'', ''A3'') : ';
g = om.get_soln('lin', 'l_Ax', 'A3');
t_is(g, f(3), 14, [t 'l3 - A3 * x']);

t = 'parse_soln : ';
t_is(om.soln.var.val.x, om.get_soln('var', 'x'), 14, [t 'var.val.x']);

t_end;
