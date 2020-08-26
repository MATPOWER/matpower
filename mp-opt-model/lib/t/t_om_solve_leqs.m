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

n = 6;

t_begin(n*length(cfg), quiet);

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
    
    t = sprintf('%s - sparse A matrix : ', name);
    om = opt_model;
    om.add_var('x', 3, zeros(size(x2)));
    om.add_lin_constraint('A12', A2(1:2, :), b2(1:2), b2(1:2));
    om.add_lin_constraint('A3', A2(3, :), b2(3), b2(3));
    opt = struct('leq_opt', struct('solver', alg));
    [x, f, e, out, J] = om.solve(opt);
    t_is(x, x2, 14, [t 'x']);
    t_is(f, A2*x-b2, 14, [t 'f']);
    t_is(e, 1, 14, [t 'exitflag']);
    t_ok(strcmp(out.alg, alg), [t 'output']);
    t_is(J, A2, 14, [t 'J']);
end

t_end;
