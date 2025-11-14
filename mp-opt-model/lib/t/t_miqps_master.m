function t_miqps_master(quiet)
% t_miqps_master - Tests of MILP/MIQP solvers via miqps_master.

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

n = 83;
nqp = 28;
nmiqp = 11;
t_begin(n*length(algs), quiet);

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
        [x, f, s, out, lam] = miqps_master([], c, A, l, u, xmin, [], [], [], opt);
        t_is(s, 1, 12, [t 'success']);
        t_is(x, [0; 15; 3], 6, [t 'x']);
        t_is(f, -78, 6, [t 'f']);
        t_is(lam.mu_l, [0;1.5;0], 9, [t 'lam.mu_l']);
        t_is(lam.mu_u, [0;0;0.5], 9, [t 'lam.mu_u']);
        t_is(lam.lower, [1;0;0], 9, [t 'lam.lower']);
        t_is(lam.upper, zeros(size(x)), 9, [t 'lam.upper']);

        if does_qp(k)
            t = sprintf('%s - unconstrained 3-d quadratic : ', names{k});
            %% from http://www.akiti.ca/QuadProgEx0Constr.html
            H = [5 -2 -1; -2 4 3; -1 3 5];
            c = [2; -35; -47];
            x0 = [0; 0; 0];
            [x, f, s, out, lam] = miqps_master(H, c, [], [], [], [], [], [], [], opt);
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
            [x, f, s, out, lam] = miqps_master(H, c, A, l, u, xmin, [], x0, [], opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [2; 4]/3, 7, [t 'x']);
            t_is(f, -74/9, 6, [t 'f']);
            t_is(lam.mu_l, [0;0;0], 13, [t 'lam.mu_l']);
            t_is(lam.mu_u, [28;4;0]/9, 4, [t 'lam.mu_u']);
            t_is(lam.lower, zeros(size(x)), 7, [t 'lam.lower']);
            t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

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
            [x, f, s, out, lam] = miqps_master(H, c, A, l, u, xmin, [], x0, [], opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [0; 2.8; 0.2; 0]/3, 5, [t 'x']);
            t_is(f, 3.29/3, 6, [t 'f']);
            t_is(lam.mu_l, [6.58;0]/3, 6, [t 'lam.mu_l']);
            t_is(lam.mu_u, [0;0], 13, [t 'lam.mu_u']);
            t_is(lam.lower, [2.24;0;0;1.7667], 4, [t 'lam.lower']);
            t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

            t = sprintf('%s - (struct) constrained 4-d QP : ', names{k});
            p = struct('H', H, 'A', A, 'l', l, 'u', u, 'xmin', xmin, 'x0', x0, 'opt', opt);
            [x, f, s, out, lam] = miqps_master(p);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [0; 2.8; 0.2; 0]/3, 5, [t 'x']);
            t_is(f, 3.29/3, 6, [t 'f']);
            t_is(lam.mu_l, [6.58;0]/3, 6, [t 'lam.mu_l']);
            t_is(lam.mu_u, [0;0], 13, [t 'lam.mu_u']);
            t_is(lam.lower, [2.24;0;0;1.7667], 4, [t 'lam.lower']);
            t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);
        else
            t_skip(nqp, sprintf('%s does not handle QP problems', names{k}));
        end

        t = sprintf('%s - infeasible LP : ', names{k});
        p = struct('A', sparse([1 1]), 'c', [1;1], 'u', -1, 'xmin', [0;0], 'opt', opt);
        [x, f, s, out, lam] = miqps_master(p);
        t_ok(s <= 0, [t 'no success']);

% opt.verbose = 3;
        t = sprintf('%s - 2-d ILP : ', names{k});
        %% from MOSEK 6.0 Guided Tour, section  7.13.1, https://docs.mosek.com/6.0/toolbox/node009.html
        c = [-2; -3];
        A = sparse([195 273; 4 40]);
        u = [1365; 140];
        xmax = [4; Inf];
        vtype = 'I';
        p = struct('c', c, 'A', A, 'u', u, 'xmax', xmax, 'vtype', vtype, 'opt', opt);
        [x, f, s, out, lam] = miqps_master(p);
        t_is(s, 1, 12, [t 'success']);
        t_is(x, [4; 2], 12, [t 'x']);
        t_is(f, -14, 12, [t 'f']);

        t = sprintf('%s - 2-d ILP (integer relaxed) : ', names{k});
        p.opt = opt_r;
        [x, f, s, out, lam] = miqps_master(p);
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
        p = struct('c', c, 'A', A, 'l', l, 'xmin', xmin, 'xmax', xmax, 'vtype', vtype, 'opt', opt);
        [x, f, s, out, lam] = miqps_master(p);
        t_is(s, 1, 12, [t 'success']);
        t_ok(norm(x - [1; 0; 3; 0; 0; 2], Inf) < 1e-12 || ...
             norm(x - [0; 0; 3; 1; 0; 2], Inf) < 1e-12 || ...
             norm(x - [0; 0; 3; 0; 2; 1], Inf) < 1e-12, [t 'x']);
        t_is(f, 16, 12, [t 'f']);

        if exist('mp.opt_model')
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
            mm = mp.opt_model();
            mm.var.add('u', 2, [], 0, 1, 'B');
            mm.var.init_indexed_name('y', {2});
            mm.var.init_indexed_name('z', {2});
            for p = 1:2
                mm.var.add('y', {p}, 3, [], 0);     % Product Y, Plant p to Customers 1-3
                mm.var.add('z', {p}, 3, [], 0);     % Product Z, Plant p to Customers 1-3
            end
            mm.qdc.add(mm.var, 'fixed', [], PlantFixedCost(:, Scenario), [], {'u'});
            mm.qdc.init_indexed_name('delivery_y', {2});
            mm.qdc.init_indexed_name('delivery_z', {2});
            for p = 1:2
                vs_y = struct('name', 'y', 'idx', {{p}});    % variable set for 'y{p}'
                vs_z = struct('name', 'z', 'idx', {{p}});    % variable set for 'z{p}'
                mm.qdc.add(mm.var, 'delivery_y', {p}, [], DeliveryCostY(:, p), [], vs_y);
                mm.qdc.add(mm.var, 'delivery_z', {p}, [], DeliveryCostZ(:, p), [], vs_z);
            end
            Au = -spdiags(PlantCapacity, 0, 2, 2);
            Ayz = sparse([1 1 1 1 1 1 0 0 0 0 0 0;   % sum for plant 1
                        0 0 0 0 0 0 1 1 1 1 1 1]); % sum for plant 2
            ub = 0;    % constraint upper bound (no lower bound)
            mm.lin.add(mm.var, 'capacity', [Au Ayz], [], ub);
            vs_y = struct('name', 'y', 'idx', {{1}, {2}});
            Ad = [speye(3) speye(3)];
            mm.lin.add(mm.var, 'demand_y', Ad, CustomerDemandY, CustomerDemandY, vs_y);
            vs_z = struct('name', 'z', 'idx', {{1}, {2}});
            mm.lin.add(mm.var, 'demand_z', Ad, CustomerDemandZ, CustomerDemandZ, vs_z);
            ef = [  490 1130/3 540;
                    410 1000/3 440;
                    410 317 410 ];
            ex = [  1 0 10 15 5 8 12 6 0 0 0 0 0 0;
                    0 1 0 0 0 0 0 0 10 15 5 8 12 6;
                    1 1 10 15 5 0 0 0 0 0 0 8 12 6;
                    0.5 1.3/3 10 15 5 0 0 0 0 0 0 8 12 6]';
            for Scenario = 1:3
                t = sprintf('%s - 14-d MILP Scenario %d: ', names{k}, Scenario);
                mm.qdc.set_params(mm.var, 'fixed', 'c', PlantFixedCost(:, Scenario));
                [x0, xl, xu, vt] = mm.var.params();
                [A, l, u] = mm.lin.params(mm.var);
                [HH, cc, kk] = mm.qdc.params(mm.var);

                [x, f, s, out, lam] = miqps_master(HH, cc, A, l, u, xl, xu, x0, vt, opt);
                t_is(s, 1, 12, [t 'exitflag']);
                t_is(x, ex(:, Scenario), 12, [t 'x']);
                t_is(f, ef(Scenario, 1), 12, [t 'f']);

                t = sprintf('%s - 14-d MILP Scenario %d (integer relaxed) : ', names{k}, Scenario);
                [x, f, s, out, lam] = miqps_master(HH, cc, A, l, u, xl, xu, x0, vt, opt_r);
                t_is(s, 1, 12, [t 'exitflag']);
                t_is(x, ex(:, 4), 12, [t 'x']);
                t_is(f, ef(Scenario, 2), 12, [t 'f']);

                t = sprintf('%s - 14-d MILP Scenario %d (integer fixed) : ', names{k}, Scenario);
                x0 = ones(size(x));
                [x, f, s, out, lam] = miqps_master(HH, cc, A, l, u, xl, xu, x0, vt, opt_f);
                t_is(s, 1, 12, [t 'exitflag']);
                t_is(x, ex(:, 3), 12, [t 'x']);
                t_is(f, ef(Scenario, 3), 12, [t 'f']);
            end
        else
            t_skip(27, sprintf('%s : mp.opt_model not available', names{k}));
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
            p = struct('H', H, 'c', c, 'A', A, 'l', l, 'u', u, ...
                'xmin', xmin, 'xmax', xmax, 'vtype', vtype, 'opt', opt);
            [x, f, s, out, lam] = miqps_master(p);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [7; 7; 0; 2], 7, [t 'x']);
            t_is(f, 1618.5, 4, [t 'f']);
            t_is(lam.mu_l, [466; 0; 0], 6, [t 'lam.mu_l']);
            t_is(lam.mu_u, [0; 272; 0], 6, [t 'lam.mu_u']);
            t_is(lam.lower, [0; 0; 349.5; 4350], 5, [t 'lam.lower']);
            t_is(lam.upper, [0; 0; 0; 0], 7, [t 'lam.upper']);

            t = sprintf('%s - 6-d IQP : ', names{k});
            %% from Bragin, et. al. https://doi.org/10.1007/s10957-014-0561-3
            %% with sign of A(2,[2;4;6]) corrected
            H = sparse(1:6, 1:6, [1 0.2 1 0.2 1 0.2], 6, 6);
            A = [-1 0.2 -1 0.2 -1 0.2;
                 -5  -1 -5  -1 -5  -1];
            u = [-48; -250];
            xmin = zeros(6, 1);
            vtype = 'I';
            p = struct('H', H, 'A', A, 'u', u, ...
                'xmin', xmin, 'vtype', vtype, 'opt', opt);
            [x, f, s, out, lam] = miqps_master(p);
            t_is(s, 1, 12, [t 'success']);
            t_ok(norm(x([1;3;5]) - [16; 16; 17], Inf) < 1e-5 || ...
                 norm(x([1;3;5]) - [17; 16; 16], Inf) < 1e-5, [t 'x([1;3;5])']);
            t_ok(norm(x([2;4;6]) - [1; 2; 2], Inf) < 1e-5 || ...
                 norm(x([2;4;6]) - [2; 2; 1], Inf) < 1e-5, [t 'x([2;4;6])']);
            t_is(f, 401.4, 5, [t 'f']);
        else
            t_skip(nmiqp, sprintf('%s does not handle MIQP problems', names{k}));
        end
% opt.verbose = 0;
    end
end

if have_feature('quadprog') && have_feature('quadprog', 'vnum') == 7.005
    warning(s1.state, diff_alg_warn_id);
end

t_end;

function TorF = have_milp_solver()
TorF = have_feature('cplex') || have_feature('highs') || ...
    have_feature('glpk') || have_feature('gurobi') || ...
    have_feature('intlinprog') || have_feature('mosek');
