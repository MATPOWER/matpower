function t_mm_solve_qps(quiet)
% t_mm_solve_qps - Tests of QP solvers via mp.opt_model.solve.

%   MP-Opt-Model
%   Copyright (c) 2010-2025, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin < 1
    quiet = 0;
end

algs = {'DEFAULT', 'BPMPD', 'MIPS', 250, 'IPOPT', 'OT', 'HIGHS', 'CPLEX', 'MOSEK', 'GUROBI', 'CLP', 'GLPK', 'OSQP', 'KNITRO'};
names = {'DEFAULT', 'BPMPD_MEX', 'MIPS', 'sc-MIPS', 'IPOPT', 'linprog/quadprog', 'HiGHS', 'CPLEX', 'MOSEK', 'Gurobi', 'CLP', 'glpk', 'OSQP', 'Knitro'};
check = {[], 'bpmpd', [], [], 'ipopt', 'quadprog', 'highs', 'cplex', 'mosek', 'gurobi', 'clp', 'glpk', 'osqp', 'knitro'};
does_qp = [1 1 1 1 1 1 1 1 1 1 1 0 1 1];

n = 30;
nqp = 21;
diff_tool = 'bbdiff';
show_diff_on_fail = false;

reps = {};

t_begin(29+n*length(algs), quiet);

diff_alg_warn_id = 'optim:linprog:WillRunDiffAlg';
if have_feature('quadprog') && have_feature('quadprog', 'vnum') == 7.005
    s1 = warning('query', diff_alg_warn_id);
    warning('off', diff_alg_warn_id);
end

for k = 1:length(algs)
    if ~isempty(check{k}) && ~have_feature(check{k})
        t_skip(n, sprintf('%s not installed', names{k}));
    else
        opt = struct('verbose', 0, 'alg', algs{k});
        mpopt = struct( ...
            'verbose', 0, ...
            'opf', struct( ...
                'violation', 1e-8 ), ...
            'cplex', struct( ...
                'lpmethod', 0, ...
                'qpmethod', 0, ...
                'opt', 0 ), ...
            'mosek', struct( ...
                'lp_alg', 0, ...
                'max_it', 0, ...
                'gap_tol', 0, ...
                'max_time', 0, ...
                'num_threads', 0, ...
                'opt', 0 ), ...
            'knitro', struct( ...
                 'tol_x', 0, ...
                 'tol_f', 0, ...
                 'maxit', 0, ...
                 'opt', 0) ...
        );
        opt.mips_opt.comptol = 1e-8;
%         if have_feature('linprog')
%             opt.linprog_opt.Algorithm = 'interior-point';
%             opt.linprog_opt.Algorithm = 'active-set';
%             opt.linprog_opt.Algorithm = 'simplex';
%             opt.linprog_opt.Algorithm = 'dual-simplex';
%         end
%         if have_feature('highs')
%             opt.highs_opt.primal_feasibility_tolerance = 1e-10;
%             opt.highs_opt.dual_feasibility_tolerance = 1e-10;
%             opt.highs_opt.ipm_optimality_tolerance = 1e-12;
%             opt.highs_opt.primal_residual_tolerance = 1e-10;
%             opt.highs_opt.dual_residual_tolerance = 1e-10;
%         end
        if have_feature('cplex')
            % alg = 0;        %% default uses barrier method with NaN bug in lower lim multipliers
            alg = 2;        %% use dual simplex
            mpopt.cplex.lpmethod = alg;
            mpopt.cplex.qpmethod = min([4 alg]);
            opt.cplex_opt = cplex_options([], mpopt);
        end
        if have_feature('mosek')
%             sc = mosek_symbcon;
%             alg = sc.MSK_OPTIMIZER_DUAL_SIMPLEX;    %% use dual simplex
%             alg = sc.MSK_OPTIMIZER_INTPNT;          %% use interior point
%             mpopt.mosek.lp_alg = alg;
            mpopt.mosek.gap_tol = 1e-10;
%             mpopt.mosek.opts.MSK_DPAR_INTPNT_TOL_PFEAS = 1e-10;
%             mpopt.mosek.opts.MSK_DPAR_INTPNT_TOL_DFEAS = 1e-10;
%             mpopt.mosek.opts.MSK_DPAR_INTPNT_TOL_INFEAS = 1e-10;
%             mpopt.mosek.opts.MSK_DPAR_INTPNT_TOL_REL_GAP = 1e-10;
            vnum = have_feature('mosek', 'vnum');
            if vnum >= 8
%                 mpopt.mosek.opts.MSK_DPAR_INTPNT_QO_TOL_PFEAS = 1e-10;
%                 mpopt.mosek.opts.MSK_DPAR_INTPNT_QO_TOL_DFEAS = 1e-10;
%                 mpopt.mosek.opts.MSK_DPAR_INTPNT_QO_TOL_INFEAS = 1e-10;
%                 mpopt.mosek.opts.MSK_DPAR_INTPNT_QO_TOL_MU_RED = 1e-10;
                mpopt.mosek.opts.MSK_DPAR_INTPNT_QO_TOL_REL_GAP = 1e-10;
            end
            opt.mosek_opt = mosek_options([], mpopt);
        end
        if have_feature('osqp')
            opt.osqp_opt.polish = 1;
%             opt.osqp_opt.alpha = 1;
%             opt.osqp_opt.eps_abs = 1e-8;
%             opt.osqp_opt.eps_rel = 1e-10;
%             opt.osqp_opt.eps_prim_inf = 1e-8;
%             opt.osqp_opt.eps_dual_inf = 1e-8;
        end
        if have_feature('knitro')
            opt.knitro_opt = artelys_knitro_options([],  mpopt);
            opt.knitro_opt.opttol = 1e-8;
        end

        t = sprintf('%s - 3-d LP : ', names{k});
        %% based on example from 'doc linprog'
        c = [-5; -4; -6];
        A = [1 -1  1;
             -3  -2  -4;
             3  2  0];
        l = [-Inf; -42; -Inf];
        u = [20; Inf; 30];
        xmin = [0; 0; 0];
        x0 = [];
        mm = mp.opt_model;
        mm.var.add('x', 3, x0, xmin);
        mm.qdc.add(mm.var, 'c', [], c);
        mm.lin.add(mm.var, 'Ax', A, l, u);
        [x, f, s, out, lam] = mm.solve(opt);
        t_is(s, 1, 12, [t 'success']);
        t_is(x, [0; 15; 3], 6, [t 'x']);
        t_is(f, -78, 6, [t 'f']);
        t_is(lam.mu_l, [0;1.5;0], 9, [t 'lam.mu_l']);
        t_is(lam.mu_u, [0;0;0.5], 9, [t 'lam.mu_u']);
        if strcmp(algs{k}, 'CLP') && ~have_feature('opti_clp')
            t_skip(2, [t 'lam.lower/upper : MEXCLP does not return multipliers on var bounds']);
        else
            t_is(lam.lower, [1;0;0], 9, [t 'lam.lower']);
            t_is(lam.upper, zeros(size(x)), 9, [t 'lam.upper']);
        end
        t_ok(~mm.has_parsed_soln(), [t 'has_parsed_soln() is false']);

        if does_qp(k)
            t = sprintf('%s - unconstrained 3-d quadratic : ', names{k});
            %% from http://www.akiti.ca/QuadProgEx0Constr.html
            H = [5 -2 -1; -2 4 3; -1 3 5];
            c = [2; -35; -47];
            x0 = [0; 0; 0];
            mm = mp.opt_model;
            mm.var.add('x', 3, x0);
            mm.qdc.add(mm.var, 'c', H, c);
            [x, f, s, out, lam] = mm.solve(opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [3; 5; 7], 6.5, [t 'x']);
            t_is(f, -249, 11, [t 'f']);
            t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
            t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
            t_is(lam.lower, zeros(size(x)), 13, [t 'lam.lower']);
            t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

            t = sprintf('%s - constrained 2-d QP : ', names{k});
            %% example from 'doc quadprog'
            H = [   1   -1;
                    -1  2   ];
            c = [-2; -6];
            A = [   1   1;
                    -1  2;
                    2   1   ];
            l = [];
            u = [2; 2; 3];
            xmin = [0; 0];
            x0 = [];
            mm = mp.opt_model;
            mm.var.add('x', 2, x0, xmin);
            mm.qdc.add(mm.var, 'c', H, c);
            mm.lin.add(mm.var, 'Ax', A, l, u);
            [x, f, s, out, lam] = mm.solve(opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [2; 4]/3, 7, [t 'x']);
            t_is(f, -74/9, 6, [t 'f']);
            t_is(lam.mu_l, [0;0;0], 13, [t 'lam.mu_l']);
            t_is(lam.mu_u, [28;4;0]/9, 4, [t 'lam.mu_u']);
            if strcmp(algs{k}, 'CLP') && ~have_feature('opti_clp')
                t_skip(2, [t 'lam.lower/upper : MEXCLP does not return multipliers on var bounds']);
            else
                t_is(lam.lower, zeros(size(x)), 7, [t 'lam.lower']);
                t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);
            end

            t = sprintf('%s - constrained 4-d QP : ', names{k});
            %% from https://v8doc.sas.com/sashtml/iml/chap8/sect12.htm
            H = [   1003.1  4.3     6.3     5.9;
                    4.3     2.2     2.1     3.9;
                    6.3     2.1     3.5     4.8;
                    5.9     3.9     4.8     10  ];
            c = zeros(4,1);
            A = [   1       1       1       1;
                    0.17    0.11    0.10    0.18    ];
            l = [1; 0.10];
            u = [1; Inf];
            xmin = zeros(4,1);
            x0 = [1; 0; 0; 1];
            mm = mp.opt_model;
            mm.var.add('x', 4, x0, xmin);
            mm.qdc.add(mm.var, 'c', H, c);
            mm.lin.add(mm.var, 'Ax', A, l, u);
            [x, f, s, out, lam] = mm.solve(opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [0; 2.8; 0.2; 0]/3, 5, [t 'x']);
            t_is(f, 3.29/3, 6, [t 'f']);
            t_is(lam.mu_l, [6.58;0]/3, 6, [t 'lam.mu_l']);
            t_is(lam.mu_u, [0;0], 13, [t 'lam.mu_u']);
            if strcmp(algs{k}, 'CLP') && ~have_feature('opti_clp')
                t_skip(2, [t 'lam.lower/upper : MEXCLP does not return multipliers on var bounds']);
            else
                t_is(lam.lower, [2.24;0;0;1.7667], 4, [t 'lam.lower']);
                t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);
            end
        else
            t_skip(nqp, sprintf('%s does not handle QP problems', names{k}));
        end

        t = sprintf('%s - infeasible LP : ', names{k});
        p = struct('A', sparse([1 1]), 'c', [1;1], 'u', -1, 'xmin', [0;0], 'opt', opt);
        [x, f, s, out, lam] = qps_master(p);
        t_ok(s <= 0, [t 'no success']);
    end
end

t = 'mm.soln.';
opt.alg = 'MIPS';
%% from https://v8doc.sas.com/sashtml/iml/chap8/sect12.htm
H = [   1003.1  4.3     6.3     5.9;
        4.3     2.2     2.1     3.9;
        6.3     2.1     3.5     4.8;
        5.9     3.9     4.8     10  ];
c = zeros(4,1);
A = [   1       1       1       1;
        0.17    0.11    0.10    0.18    ];
l = [1; 0.10];
u = [1; Inf];
xmin = zeros(4,1);
x0 = [1; 0; 0; 1];
mm = mp.opt_model;
mm.var.add('x', 4, x0, xmin);
mm.qdc.add(mm.var, 'c', H, c);
mm.lin.add(mm.var, 'Ax', A, l, u);
opt.parse_soln = 1;
[x, f, s, out, lam] = mm.solve(opt);
t_is(mm.soln.x, x, 14, [t 'x']);
t_is(mm.soln.f, f, 14, [t 'f']);
t_is(mm.soln.eflag, s, 14, [t 'eflag']);
t_str_match(mm.soln.output.alg, out.alg, [t 'output.alg']);
t_is(mm.soln.lambda.lower, lam.lower, 14, [t 'mm.soln.lambda.lower']);
t_is(mm.soln.lambda.upper, lam.upper, 14, [t 'mm.soln.lambda.upper']);
t_is(mm.soln.lambda.mu_l, lam.mu_l, 14, [t 'mm.soln.lambda.mu_l']);
t_is(mm.soln.lambda.mu_u, lam.mu_u, 14, [t 'mm.soln.lambda.mu_u']);

t = 'mm.var.get_soln(mm.soln, ''x'') : ';
[x1, mu_l, mu_u] = mm.var.get_soln(mm.soln, 'x');
t_is(x1, x, 14, [t 'x']);
t_is(mu_l, lam.lower, 14, [t 'mu_l']);
t_is(mu_u, lam.upper, 14, [t 'mu_u']);

t = 'mm.var.get_soln(mm.soln, ''mu_l'', ''x'') : ';
t_is(mm.var.get_soln(mm.soln, 'mu_l', 'x'), lam.lower, 14, [t 'mu_l']);

t = 'mm.lin.get_soln(mm.var, mm.soln, ''Ax'') : ';
[g, mu_l] = mm.lin.get_soln(mm.var, mm.soln, 'Ax');
t_is(g{1}, A*x-u, 14, [t 'A * x - u']);
t_is(g{2}, l-A*x, 14, [t 'l - A * x']);
t_is(mu_l, lam.mu_l, 14, [t 'mu_l']);

t = 'mm.lin.get_soln(mm.var, mm.soln, {''mu_u'', ''mu_l'', ''Ax_u''}, ''Ax'') : ';
[mu_u, mu_l, g] = mm.lin.get_soln(mm.var, mm.soln, {'mu_u', 'mu_l', 'Ax_u'}, 'Ax');
t_is(g, A*x-u, 14, [t 'A * x - u']);
t_is(mu_l, lam.mu_l, 14, [t 'mu_l']);
t_is(mu_u, lam.mu_u, 14, [t 'mu_u']);

t = 'mm.qdc.get_soln(mm.var, mm.soln, ''c'') : ';
[f1, df, d2f] = mm.qdc.get_soln(mm.var, mm.soln, 'c');
t_is(f1, f, 14, [t 'f']);
t_is(df, H*x+c, 14, [t 'df']);
t_is(d2f, H, 14, [t 'd2f']);

t = 'mm.qdc.get_soln(mm.var, mm.soln, ''df'', ''c'') : ';
df = mm.qdc.get_soln(mm.var, mm.soln, 'df', 'c');
t_is(df, H*x+c, 14, [t 'df']);

t = 'parse_soln : ';
t_ok(mm.has_parsed_soln(), [t 'has_parsed_soln() is true']);
t_is(mm.var.soln.val.x, mm.var.get_soln(mm.soln, 'x'), 14, [t 'var.val.x']);
t_is(mm.var.soln.mu_l.x, mm.var.get_soln(mm.soln, 'mu_l', 'x'), 14, [t 'var.mu_l.x']);
t_is(mm.var.soln.mu_u.x, mm.var.get_soln(mm.soln, 'mu_u', 'x'), 14, [t 'var.mu_u.x']);
t_is(mm.lin.soln.mu_l.Ax, mm.lin.get_soln(mm.var, mm.soln, 'mu_l', 'Ax'), 14, [t 'lin.mu_l.Ax']);
t_is(mm.lin.soln.mu_u.Ax, mm.lin.get_soln(mm.var, mm.soln, 'mu_u', 'Ax'), 14, [t 'lin.mu_u.Ax']);

t = 'disp_soln';
if have_feature('matlab') && have_feature('matlab', 'vnum') < 8
    t_skip(1, [t ' : old MATLAB number format not exact match']);
else
    rn = fix(1e9*rand);
    [pathstr, name, ext] = fileparts(which('t_opt_model'));
    fname = 't_mm_solve_qps_display_soln';
    fname_e = fullfile(pathstr, 'display_soln', sprintf('%s.txt', fname));
    fname_g = sprintf('%s_%d.txt', fname, rn);
    [fd, msg] = fopen(fname_g, 'wt');   %% open solution file
    if fd == -1
        error('t_mm_solve_qps: could not create %d : %s', fname, msg);
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
end

if have_feature('quadprog') && have_feature('quadprog', 'vnum') == 7.005
    warning(s1.state, diff_alg_warn_id);
end

t_end;
