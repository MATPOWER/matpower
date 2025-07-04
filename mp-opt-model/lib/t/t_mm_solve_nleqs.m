function t_mm_solve_nleqs(quiet)
% t_mm_solve_nleqs - Tests of NLEQ solvers via mp.opt_model.solve.

%   MP-Opt-Model
%   Copyright (c) 2010-2025, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin < 1
    quiet = 0;
end

core_sp_newton = struct( ...
    'alg',              'NEWTON', ...
    'name',             'Newton''s', ...
    'default_max_it',   10, ...
    'need_jac',         1, ...
    'update_fcn',       @(x, f, J)newton_update_fcn(x, f, J, '')  );
core_sp_gs = struct( ...
    'alg',              'GS', ...
    'name',             'Gauss-Seidel', ...
    'default_max_it',   1000, ...
    'need_jac',         0, ...
    'update_fcn',       @(x, f)x_update_fcn2(x, f)  );

if have_feature('matlab')
    %%  alg         name        check       opts
    cfg = {
        {'DEFAULT', 'default',  []          []  },
        {'NEWTON',  'Newton',   [],         []  },
        {'FD',      'fast-decoupled Newton',[],[]  },
        {'FSOLVE',  'fsolve-1', 'fsolve',   []  },
        {'FSOLVE',  'fsolve-2', 'fsolve',   struct('Algorithm', 'trust-region-dogleg')  },
        {'FSOLVE',  'fsolve-3', 'fsolve',   struct('Algorithm', 'trust-region-reflective')  },
        {'FSOLVE',  'fsolve-4', 'fsolve',   struct('Algorithm', 'levenberg-marquardt', 'TolX', 1e-11) },
        {'GS',      'Gauss-Seidel',[],      []  },
        {'CORE-N',  'Newton-CORE',   [],    core_sp_newton  },
        {'CORE-GS', 'Gauss-Seidel-CORE',[], core_sp_gs  },
    };
    if have_feature('matlab', 'vnum') <= 7.010
        cfg([6]) = [];  %% MATLAB 7.10 does not work w/ fsolve alg 3
    end
else    %% octave
    %%  alg         name        check       opts
    cfg = {
        {'DEFAULT', 'default',  []          []  },
        {'NEWTON',  'Newton',   [],         []  },
        {'FD',      'fast-decoupled Newton',[],[]  },
        {'FSOLVE',  'fsolve', 'fsolve',     struct('TolX', 1e-11)  },
        {'GS',      'Gauss-Seidel',[],      []  },
        {'CORE-N',  'Newton-CORE',   [],    core_sp_newton  },
        {'CORE-GS', 'Gauss-Seidel-CORE',[], core_sp_gs  },
    };
end

n = 18;
diff_tool = 'bbdiff';
show_diff_on_fail = false;

reps = {};

t_begin(16+n*length(cfg), quiet);

for k = 1:length(cfg)
    alg   = cfg{k}{1};
    name  = cfg{k}{2};
    check = cfg{k}{3};
    opts  = cfg{k}{4};
    if ~isempty(check) && ~have_feature(check)
        t_skip(n, sprintf('%s not installed', name));
    else
        opt = struct('verbose', 0, 'alg', alg, 'tol', 1e-11);
        switch alg
            case {'DEFAULT', 'NEWTON'}
            case {'FD'}
                opt.fd_opt.jac_approx_fcn = @jac_approx_fcn2;
            case 'FSOLVE'
                opt.fsolve_opt = opts;
            case 'GS'
                opt.gs_opt.x_update_fcn = @(x, f)x_update_fcn2(x, f);
            case {'CORE-N', 'CORE-GS'}
                opt.core_sp = opts;
                opt.alg = 'CORE';
        end

        switch alg
            case {'DEFAULT', 'NEWTON', 'FSOLVE', 'CORE-N'}
                t = sprintf('%s - 2-d function : ', name);
                x0 = [-1;0];
                mm = mp.opt_model;
                mm.var.add('x', 2, x0);
                mm.nle.add(mm.var, 'f', 2, @f1, []);
                [x, f, e, out, jac] = mm.solve(opt);
                t_is(e, 1, 12, [t 'success']);
                t_is(x, [-3; 4], 8, [t 'x']);
                t_is(f, 0, 10, [t 'f']);
                switch alg
                    case {'DEFAULT', 'CORE-N'}
                        out_alg = 'NEWTON';
                    otherwise
                        out_alg = alg;
                end
                t_str_match(out.alg, out_alg, [t 'out.alg']);
                eJ = [1 1; 6 1];
                t_is(jac, eJ, 5.8, [t 'jac']);

                t = sprintf('%s - 2-d function (max_it) : ', name);
                opt.max_it = 3;
                [x, f, e, out] = mm.solve(opt);
                t_is(e, 0, 12, [t 'no success']);
                t_ok(out.iterations == 3 || out.iterations == 4, [t 'iterations']);
                opt = rmfield(opt, 'max_it');

                t = sprintf('%s - 5-d lin/nonlin function : ', name);
                x0_1 = [-1;0];
                x0_2 = [0;0;0];
                A2 = sparse([2 -1 0; -3 1 -2; 0 5 -4]);
                b2 = [-5; 1; -7];
                x2 = [-2; 1; 3];
                mm = mp.opt_model;
                mm.var.add('x1', 2, x0_1);
                mm.var.add('x2', 3, x0_2);
                mm.nle.add(mm.var, 'f', 2, @f1, [], {'x1'});
                mm.lin.add(mm.var, 'Ax_b', A2, b2, b2, {'x2'});
                [x, f, e, out, jac] = mm.solve(opt);
                t_is(e, 1, 12, [t 'success']);
                t_is(x, [-3; 4; -2; 1; 3], 8, [t 'x']);
                t_is(f, 0, 10, [t 'f']);
                t_str_match(out.alg, out_alg, [t 'out.alg']);
                eJ = [[1 1; 6 1] zeros(2, 3);
                      zeros(3, 2) A2 ];
                t_is(jac, eJ, 5.8, [t 'jac']);
            otherwise
                t_skip(12, sprintf('not implemented for solver ''%s''', alg));
        end

        t = sprintf('%s - 2-d function2 (struct) : ', name);
        x0 = [1;2];
        mm = mp.opt_model;
        mm.var.add('x', 2, x0);
        mm.nle.add(mm.var, 'f', 2, @f2, []);
        [x, f, e] = mm.solve(opt);
        t_is(e, 1, 12, [t 'success']);
        t_is(x, [2; 3], 8, [t 'x']);
        t_is(f, 0, 10, [t 'f']);
        t_ok(~mm.has_parsed_soln(), [t 'has_parsed_soln() is false']);

        opt.max_it = 3;
        t = sprintf('%s - 2-d function2 (max_it) : ', name);
        [x, f, e, out] = mm.solve(opt);
        t_is(e, 0, 12, [t 'no success']);
        t_ok(out.iterations == 3 || out.iterations == 4, [t 'iterations']);
        opt = rmfield(opt, 'max_it');
    end
end

t = 'mm.soln.';
opt.alg = 'DEFAULT';
x0_1 = [-1;0];
x0_2 = [0;0;0];
A2 = sparse([2 -1 0; -3 1 -2; 0 5 -4]);
b2 = [-5; 1; -7];
x2 = [-2; 1; 3];
mm = mp.opt_model;
mm.var.add('x1', 2, x0_1);
mm.var.add('x2', 3, x0_2);
mm.nle.add(mm.var, 'f', 2, @f1, [], {'x1'});
mm.lin.add(mm.var, 'Ax_b', A2, b2, b2, {'x2'});
opt.parse_soln = 1;
[x, f, e, out, jac] = mm.solve(opt);
t_is(mm.soln.x, x, 14, [t 'x']);
t_is(mm.soln.f, f, 14, [t 'f']);
t_is(mm.soln.eflag, e, 14, [t 'eflag']);
t_str_match(mm.soln.output.alg, out.alg, [t 'output.alg']);
t_is(mm.soln.jac, jac, 14, [t 'jac']);

t = 'mm.var.get_soln(mm.soln, ''x1'') : ';
t_is(mm.var.get_soln(mm.soln, 'x1'), x(1:2), 14, [t 'x1']);

t = 'mm.var.get_soln(mm.soln, ''x'', ''x2'') : ';
t_is(mm.var.get_soln(mm.soln, 'x', 'x2'), x(3:5), 14, [t 'x2']);

t = 'mm.lin.get_soln(mm.var, mm.soln, ''g'', ''Ax_b'') : ';
g = mm.lin.get_soln(mm.var, mm.soln, 'g', 'Ax_b');
t_is(g{1}, f(3:5), 14, [t 'A * x - u']);
t_is(g{2}, f(3:5), 14, [t 'l - A * x']);

t = 'mm.nle.get_soln(mm.var, mm.soln, true, ''f'') : ';
g = mm.nle.get_soln(mm.var, mm.soln, true, 'f');
t_is(g, f(1:2), 14, [t 'f']);

t = 'mm.nle.get_soln(mm.var, mm.soln, true, {''g'', ''dg''}, ''f'') : ';
[g, dg] = mm.nle.get_soln(mm.var, mm.soln, true, {'g', 'dg'}, 'f');
t_is(g, f(1:2), 14, [t 'f']);
t_is(dg, jac(1:2, 1:2), 14, [t 'jac']);

t = 'parse_soln : ';
t_ok(mm.has_parsed_soln(), [t 'has_parsed_soln() is true']);
t_is(mm.var.soln.val.x1, mm.var.get_soln(mm.soln, 'x1'), 14, [t 'var.val.x1']);
t_is(mm.var.soln.val.x2, mm.var.get_soln(mm.soln, 'x2'), 14, [t 'var.val.x2']);

t = 'disp_soln';
rn = fix(1e9*rand);
[pathstr, name, ext] = fileparts(which('t_opt_model'));
fname = 't_mm_solve_nleqs_display_soln';
fname_e = fullfile(pathstr, 'display_soln', sprintf('%s.txt', fname));
fname_g = sprintf('%s_%d.txt', fname, rn);
[fd, msg] = fopen(fname_g, 'wt');   %% open solution file
if fd == -1
    error('t_mm_solve_nleqs: could not create %d : %s', fname, msg);
end
mm.display_soln(fd);    %% write out solution
fclose(fd);
if ~t_file_match(fname_g, fname_e, t, reps, 1);
    fprintf('  compare these 2 files:\n    %s\n    %s\n', fname_g, fname_e);
    if show_diff_on_fail
        cmd = sprintf('%s %s %s', diff_tool, fname_g, fname_e);
        [status, result] = system(cmd);
        keyboard
    end
end

t_end;


%% 2-d problem with 2 solutions
%% from https://www.chilimath.com/lessons/advanced-algebra/systems-non-linear-equations/
function [f, J] = f1(x)
if iscell(x)    %% in case it was part of a varset
    x = x{1};
end
f = [  x(1)   + x(2) - 1;
      -x(1)^2 + x(2) + 5    ];
if nargout > 1
    J = sparse([1 1; -2*x(1) 1]);
end

%% another 2-d problem
%% from Christi Patton Luks, https://www.youtube.com/watch?v=pJG4yhtgerg
function [f, J] = f2(x)
f = [  x(1)^2 +   x(1)*x(2)   - 10;
       x(2)   + 3*x(1)*x(2)^2 - 57  ];
if nargout > 1
    J = sparse([    2*x(1)+x(2)    x(1);
                    3*x(2)^2       6*x(1)*x(2)+1    ]);
end

function JJ = jac_approx_fcn2()
%% for use with fast-decoupled Newton's method
J = [7 2; 27 37];
JJ = {J(1,1), J(2,2)};

function x = x_update_fcn2(x, f)
%% for use with Gauss-Seidel method
x(1) = sqrt(10 - x(1)*x(2));
x(2) = sqrt((57-x(2))/3/x(1));

function x = newton_update_fcn(x, f, J, lin_solver)
dx = mplinsolve(J, -f, lin_solver);     %% compute update step
x = x + dx;                             %% update x
