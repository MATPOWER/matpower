function t_om_solve_miqps(quiet)
% t_om_solve_miqps - Tests of MIQP solvers via opt_model.solve.

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

algs = {'DEFAULT', 'CPLEX', 'MOSEK', 'GUROBI', 'HIGHS', 'GLPK', 'OT'};
names = {'DEFAULT', 'CPLEX', 'MOSEK', 'Gurobi', 'HiGHS', 'glpk', 'intlin/lin/quadprog'};
check = {@have_milp_solver, 'cplex', 'mosek', 'gurobi', 'highs', 'glpk', 'intlinprog'};
does_qp   = [0 1 1 1 1 0 1];
does_miqp = [0 1 1 1 0 0 0];
if have_feature('gurobi') || have_feature('cplex') || have_feature('mosek') ...
        || have_feature('highs') || have_feature('quadprog')
    does_qp(1) = 1;
end
if have_feature('gurobi') || have_feature('cplex') || have_feature('mosek')
    does_miqp(1) = 1;
end

n = 60;
nmiqp = 20;
nmiqp_soln = 30;
diff_tool = 'bbdiff';
show_diff_on_fail = false;

reps = {};

t_begin(nmiqp_soln+n*length(algs), quiet);

diff_alg_warn_id = 'optim:linprog:WillRunDiffAlg';
if have_feature('quadprog') && have_feature('quadprog', 'vnum') == 7.005
    s1 = warning('query', diff_alg_warn_id);
    warning('off', diff_alg_warn_id);
end

for k = 1:length(algs)
    if ~isempty(check{k}) && ...
            (ischar(check{k}) && ~have_feature(check{k}) || ...
             isa(check{k}, 'function_handle') && ~check{k}())
        t_skip(n, sprintf('%s not installed', names{k}));
    else
        opt = struct('verbose', 0, 'alg', algs{k});
        mpopt = struct( ...
            'verbose', 0, ...
            'opf', struct( ...
                'violation', 5e-6 ), ...
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
                'opt', 0 ) ...
        );
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
        opt_r = opt;
        opt_r.relax_integer = 1;
        opt_f = opt;
        opt_f.fix_integer = 1;

% opt.verbose = 3;
        t = sprintf('%s - 2-d ILP : ', names{k});
        %% from MOSEK 6.0 Guided Tour, section  7.13.1, https://docs.mosek.com/6.0/toolbox/node009.html
        c = [-2; -3];
        A = sparse([195 273; 4 40]);
        u = [1365; 140];
        xmax = [4; Inf];
        vtype = 'I';
        om = opt_model;
        om.add_var('x', 2, [], [], xmax, vtype);
        om.add_quad_cost('c', [], c);
        om.add_lin_constraint('Ax', A, [], u);
        [x, f, s, out, lam] = om.solve(opt);
        t_is(s, 1, 12, [t 'success']);
        t_is(x, [4; 2], 12, [t 'x']);
        t_is(f, -14, 12, [t 'f']);
        t_ok(~om.has_parsed_soln(), [t 'has_parsed_soln() is false']);
        t = sprintf('%s - 2-d ILP (integer relaxed) : ', names{k});
        [x, f, s, out, lam] = om.solve(opt_r);
        t_is(s, 1, 12, [t 'success']);
        t_is(x, [2.441860465; 3.255813953], 8, [t 'x']);
        t_is(f, -14.651162791, 8, [t 'f']);

        t = sprintf('%s - 6-d ILP : ', names{k});
        %% from https://doi.org/10.1109/TASE.2020.2998048
        c = [1; 2; 3; 1; 2; 3];
        A = [1 3 5 1 3 5;
             2 1.5 5 2 0.5 1];
        l = [26; 16];
        xmin = zeros(6, 1);
        xmax = 3 * ones(6, 1);
        vtype = 'I';
        om = opt_model;
        om.add_var('x', 6, [], xmin, xmax, vtype);
        om.add_quad_cost('c', [], c);
        om.add_lin_constraint('Ax', A, l, u);
        [x, f, s, out, lam] = om.solve(opt);
        t_is(s, 1, 12, [t 'success']);
        t_ok(norm(x - [1; 0; 3; 0; 0; 2], Inf) < 1e-12 || ...
             norm(x - [0; 0; 3; 1; 0; 2], Inf) < 1e-12 || ...
             norm(x - [0; 0; 3; 0; 2; 1], Inf) < 1e-12, [t 'x']);
        t_is(f, 16, 12, [t 'f']);
        t = sprintf('%s - 6-d ILP (integer relaxed) : ', names{k});
        [x, f, s, out, lam] = om.solve(opt_r);
        t_is(s, 1, 12, [t 'success']);
        t_is([x([1;2;4;5]); x(3)+x(6)], [0; 0; 0; 0; 5.2], 7, [t 'x']);
        t_is(f, 15.6, 7, [t 'f']);

        %% from milp_ex1.m
        PlantCapacity = [ 60; 60 ];     % Plants 1, 2
        CustomerDemandY = [ 10; 15; 5 ];    % Customers 1, 2, 3
        CustomerDemandZ = [  8; 12; 6 ];    % Customers 1, 2, 3
        DeliveryCostY = [
            4           7;          % Customer 1
            5           6;          % Customer 2
            3           8           % Customer 3
        ];
        DeliveryCostZ = [
            10          4;          % Customer 1
            9          5;          % Customer 2
            12          3           % Customer 3
        ];
        PlantFixedCost = [
            100         100          50;    % Plant 1
            200         100         120     % Plant 2
        ];
        Scenario = 1;   % build scenario 1 initially
        om = opt_model();
        om.add_var('u', 2, [], 0, 1, 'B');
        om.var.init_indexed_name('y', {2});
        om.var.init_indexed_name('z', {2});
        for p = 1:2
            om.var.add('y', {p}, 3, [], 0);     % Product Y, Plant p to Customers 1-3
            om.var.add('z', {p}, 3, [], 0);     % Product Z, Plant p to Customers 1-3
        end
        om.qdc.add(om.var, 'fixed', [], PlantFixedCost(:, Scenario), [], {'u'});
        om.qdc.init_indexed_name('delivery_y', {2});
        om.qdc.init_indexed_name('delivery_z', {2});
        for p = 1:2
            vs_y = struct('name', 'y', 'idx', {{p}});    % variable set for 'y{p}'
            vs_z = struct('name', 'z', 'idx', {{p}});    % variable set for 'z{p}'
            om.qdc.add(om.var, 'delivery_y', {p}, [], DeliveryCostY(:, p), [], vs_y);
            om.qdc.add(om.var, 'delivery_z', {p}, [], DeliveryCostZ(:, p), [], vs_z);
        end
        Au = -spdiags(PlantCapacity, 0, 2, 2);
        Ayz = sparse([1 1 1 1 1 1 0 0 0 0 0 0;   % sum for plant 1
                    0 0 0 0 0 0 1 1 1 1 1 1]); % sum for plant 2
        ub = 0;    % constraint upper bound (no lower bound)
        om.lin.add(om.var, 'capacity', [Au Ayz], [], ub);
        vs_y = struct('name', 'y', 'idx', {{1}, {2}});
        Ad = [speye(3) speye(3)];
        om.lin.add(om.var, 'demand_y', Ad, CustomerDemandY, CustomerDemandY, vs_y);
        vs_z = struct('name', 'z', 'idx', {{1}, {2}});
        om.lin.add(om.var, 'demand_z', Ad, CustomerDemandZ, CustomerDemandZ, vs_z);
        ef = [  490 1130/3 540;
                410 1000/3 440;
                410 317 410 ];
        ex = [  1 0 10 15 5 8 12 6 0 0 0 0 0 0;
                0 1 0 0 0 0 0 0 10 15 5 8 12 6;
                1 1 10 15 5 0 0 0 0 0 0 8 12 6;
                0.5 1.3/3 10 15 5 0 0 0 0 0 0 8 12 6]';
        for Scenario = 1:3
            t = sprintf('%s - 14-d MILP Scenario %d: ', names{k}, Scenario);
            om.qdc.set_params(om.var, 'fixed', 'c', PlantFixedCost(:, Scenario));
            [x, f, s, out, lam] = om.solve(opt);
            t_is(s, 1, 12, [t 'exitflag']);
            t_is(x, ex(:, Scenario), 12, [t 'x']);
            t_is(f, ef(Scenario, 1), 12, [t 'f']);

            t = sprintf('%s - 14-d MILP Scenario %d (integer relaxed) : ', names{k}, Scenario);
            [x, f, s, out, lam] = om.solve(opt_r);
            t_is(s, 1, 12, [t 'exitflag']);
            t_is(x, ex(:, 4), 12, [t 'x']);
            t_is(f, ef(Scenario, 2), 12, [t 'f']);
            
            t = sprintf('%s - 14-d MILP Scenario %d (integer fixed) : ', names{k}, Scenario);
            opt_f.x0 = ones(size(x));
            [x, f, s, out, lam] = om.solve(opt_f);
            t_is(s, 1, 12, [t 'exitflag']);
            t_is(x, ex(:, 3), 12, [t 'x']);
            t_is(f, ef(Scenario, 3), 12, [t 'f']);
        end

        if does_miqp(k)
            t = sprintf('%s - 4-d MIQP : ', names{k});
            %% from cplexmiqpex.m, CPLEX_Studio_Academic124/cplex/examples/src/matlab/cplexmiqpex.m
            %% Note: This is a lame example; the integer relaxed problem already
            %%       has an integer feasible solution, so this is actually just
            %%       a simple QP. -RDZ 10/29/24
            H = sparse([ 33   6    0    0;
                          6  22   11.5  0;
                          0  11.5 11    0;
                          0   0    0    0]);
            c = [-1 -2 -3 -1]';
            Aineq = [-1  1  1 10;
               1 -3  1  0];
            bineq = [20  30]';
            Aeq   = [0  1  0 -3.5];
            beq   =  0;
            xmin    = [ 0;   0;   0; 2];
            xmax    = [40; Inf; Inf; 3];
            A = sparse([Aeq; Aineq]);
            l = [beq; -Inf; -Inf];
            u = [beq; bineq];
            vtype = 'CCCI';
            om = opt_model;
            om.add_var('x', 4, [], xmin, xmax, vtype);
            om.add_quad_cost('c', H, c);
            om.add_lin_constraint('Ax', A, l, u);
            [x, f, s, out, lam] = om.solve(opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [7; 7; 0; 2], 7, [t 'x']);
            t_is(f, 1618.5, 4, [t 'f']);
            t_is(lam.mu_l, [466; 0; 0], 6, [t 'lam.mu_l']);
            t_is(lam.mu_u, [0; 272; 0], 6, [t 'lam.mu_u']);
            t_is(lam.lower, [0; 0; 349.5; 4350], 5, [t 'lam.lower']);
            t_is(lam.upper, [0; 0; 0; 0], 7, [t 'lam.upper']);
            t = sprintf('%s - 4-d MIQP (integer relaxed) : ', names{k});
            [x, f, s, out, lam] = om.solve(opt_r);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [7; 7; 0; 2], 7, [t 'x']);
            t_is(f, 1618.5, 4, [t 'f']);
            t_is(lam.mu_l, [466; 0; 0], 6, [t 'lam.mu_l']);
            t_is(lam.mu_u, [0; 272; 0], 6, [t 'lam.mu_u']);
            t_is(lam.lower, [0; 0; 349.5; 4350], 5, [t 'lam.lower']);
            t_is(lam.upper, [0; 0; 0; 0], 7, [t 'lam.upper']);

            t = sprintf('%s - 6-d IQP : ', names{k});
            %% from Bragin, et. al. https://doi.org/10.1007/s10957-014-0561-3
            H = sparse(1:6, 1:6, [1 0.2 1 0.2 1 0.2], 6, 6);
            a = [-5 1 -5 1 -5 1];
            A = [a/5; a];
            u = [-48; -250];
            xmin = zeros(6, 1);
            vtype = 'I';
            om = opt_model;
            om.add_var('x', 6, [], xmin, [], vtype);
            om.add_quad_cost('c', H, []);
            om.add_lin_constraint('Ax', A, [], u);
            [x, f, s, out, lam] = om.solve(opt);
            t_is(s, 1, 12, [t 'success']);
            t_ok(norm(x - [16; 0; 17; 0; 17; 0], Inf) < 1e-7 || ...
                 norm(x - [17; 0; 16; 0; 17; 0], Inf) < 1e-7 || ...
                 norm(x - [17; 0; 17; 0; 16; 0], Inf) < 1e-7, [t 'x']);
            t_is(f, 417, 6, [t 'f']);
            t = sprintf('%s - 6-d IQP (integer relaxed) : ', names{k});
            [x, f, s, out, lam] = om.solve(opt_r);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [50;0;50;0;50;0]/3, 8, [t 'x']);
            t_is(f, 1250/3, 6, [t 'f']);
        else
            t_skip(nmiqp, sprintf('%s does not handle MIQP problems', names{k}));
        end
% opt.verbose = 0;
    end
end

if have_milp_solver()
    t = 'om.soln.';
    c = [-2; -3];
    A = sparse([195 273; 4 40]);
    u = [1365; 140];
    xmax = [4; Inf];
    vtype = 'I';
    om = opt_model;
    om.add_var('x', 2, [], [], xmax, vtype);
    om.add_quad_cost('c', [], c);
    om.add_lin_constraint('Ax', A, [], u);
    opt.parse_soln = 1;
    [x, f, s, out, lam] = om.solve(opt);
    t_is(om.soln.x, x, 14, [t 'x']);
    t_is(om.soln.f, f, 14, [t 'f']);
    t_is(om.soln.eflag, s, 14, [t 'eflag']);
    t_str_match(om.soln.output.alg, out.alg, [t 'output.alg']);
    t_is(om.soln.lambda.lower, lam.lower, 14, [t 'om.soln.lambda.lower']);
    t_is(om.soln.lambda.upper, lam.upper, 14, [t 'om.soln.lambda.upper']);
    t_is(om.soln.lambda.mu_l, lam.mu_l, 14, [t 'om.soln.lambda.mu_l']);
    t_is(om.soln.lambda.mu_u, lam.mu_u, 14, [t 'om.soln.lambda.mu_u']);

    t = 'om.get_soln(''var'', ''x'') : ';
    [x1, mu_l, mu_u] = om.get_soln('var', 'x');
    t_is(x1, x, 14, [t 'x']);
    t_is(mu_l, lam.lower, 14, [t 'mu_l']);
    t_is(mu_u, lam.upper, 14, [t 'mu_u']);

    t = 'om.get_soln(''var'', ''mu_u'', ''x'') : ';
    t_is(om.get_soln('var', 'mu_u', 'x'), lam.upper, 14, [t 'mu_l']);

    t = 'om.get_soln(''lin'', ''Ax'') : ';
    [g, mu_u] = om.get_soln('lin', 'Ax');
    t_is(g{1}, A*x-u, 14, [t 'A * x - u']);
    t_ok(all(isinf(g{2}) & g{2} < 0), [t 'l - A * x']);
    t_is(mu_u, lam.mu_u, 14, [t 'mu_u']);

    t = 'om.get_soln(''lin'', {''mu_u'', ''mu_l'', ''g''}, ''Ax'') : ';
    [mu_u, mu_l, g] = om.get_soln('lin', {'mu_u', 'mu_l', 'g'}, 'Ax');
    t_is(g{1}, A*x-u, 14, [t 'A * x - u']);
    t_ok(all(isinf(g{2}) & g{2} < 0), [t 'l - A * x']);
    t_is(mu_l, lam.mu_l, 14, [t 'mu_l']);
    t_is(mu_u, lam.mu_u, 14, [t 'mu_u']);

    t = 'om.get_soln(''qdc'', ''c'') : ';
    [f1, df, d2f] = om.get_soln('qdc', 'c');
    t_is(sum(f1), f, 14, [t 'f']);
    t_is(df, c, 14, [t 'df']);
    t_is(d2f, zeros(2,1), 14, [t 'd2f']);

    t = 'om.get_soln(''qdc'', ''df'', ''c'') : ';
    df = om.get_soln('qdc', 'df', 'c');
    t_is(df, c, 14, [t 'df']);

    t = 'parse_soln : ';
    t_ok(om.has_parsed_soln(), [t 'has_parsed_soln() is true']);
    t_is(om.var.soln.val.x, om.get_soln('var', 'x'), 14, [t 'var.val.x']);
    t_is(om.var.soln.mu_l.x, om.get_soln('var', 'mu_l', 'x'), 14, [t 'var.mu_l.x']);
    t_is(om.var.soln.mu_u.x, om.get_soln('var', 'mu_u', 'x'), 14, [t 'var.mu_u.x']);
    t_is(om.lin.soln.mu_l.Ax, om.get_soln('lin', 'mu_l', 'Ax'), 14, [t 'lin.mu_l.Ax']);
    t_is(om.lin.soln.mu_u.Ax, om.get_soln('lin', 'mu_u', 'Ax'), 14, [t 'lin.mu_u.Ax']);

    t = 'disp_soln';
    rn = fix(1e9*rand);
    [pathstr, name, ext] = fileparts(which('t_opt_model'));
    fname = 't_mm_solve_miqps_display_soln';
    fname_e = fullfile(pathstr, 'display_soln', sprintf('%s.txt', fname));
    fname_g = sprintf('%s_%d.txt', fname, rn);
    [fd, msg] = fopen(fname_g, 'wt');   %% open solution file
    if fd == -1
        error('t_om_solve_miqps: could not create %d : %s', fname, msg);
    end
    om.display_soln(fd);    %% write out solution
    fclose(fd);
    if ~t_file_match(fname_g, fname_e, t, reps, 1);
        fprintf('  compare these 2 files:\n    %s\n    %s\n', fname_g, fname_e);
        if show_diff_on_fail
            cmd = sprintf('%s %s %s', diff_tool, fname_g, fname_e);
            [status, result] = system(cmd);
            keyboard
        end
    end
else
    t_skip(nmiqp_soln, 'no MILP/MIQP solver installed');
end

if have_feature('quadprog') && have_feature('quadprog', 'vnum') == 7.005
    warning(s1.state, diff_alg_warn_id);
end

t_end;

function TorF = have_milp_solver()
TorF = have_feature('cplex') || have_feature('highs') || ...
    have_feature('glpk') || have_feature('gurobi') || ...
    have_feature('intlinprog') || have_feature('mosek');

