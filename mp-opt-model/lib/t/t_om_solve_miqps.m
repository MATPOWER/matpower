function t_om_solve_miqps(quiet)
% t_om_solve_miqps - Tests of MIQP solvers via opt_model.solve.

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

algs = {'DEFAULT', 'CPLEX', 'MOSEK', 'GUROBI', 'GLPK', 'OT'};
names = {'DEFAULT', 'CPLEX', 'MOSEK', 'Gurobi', 'glpk', 'intlin/lin/quadprog'};
check = {@have_miqp_solver, 'cplex', 'mosek', 'gurobi', 'glpk', 'intlinprog'};
does_qp = [0 1 1 1 0 0];
if have_feature('gurobi') || have_feature('cplex') || have_feature('mosek')
    does_qp(1) = 1;
end

n = 17;
nmiqp = 10;
t_begin(29+n*length(algs), quiet);

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

        if does_qp(k)
            t = sprintf('%s - 4-d MIQP : ', names{k});
            %% from cplexmiqpex.m, CPLEX_Studio_Academic124/cplex/examples/src/matlab/cplexmiqpex.m
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
        else
            t_skip(nmiqp, sprintf('%s does not handle MIQP problems', names{k}));
        end
% opt.verbose = 0;
    end
end

if have_miqp_solver()
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
    t_is(om.soln.var.val.x, om.get_soln('var', 'x'), 14, [t 'var.val.x']);
    t_is(om.soln.var.mu_l.x, om.get_soln('var', 'mu_l', 'x'), 14, [t 'var.mu_l.x']);
    t_is(om.soln.var.mu_u.x, om.get_soln('var', 'mu_u', 'x'), 14, [t 'var.mu_u.x']);
    t_is(om.soln.lin.mu_l.Ax, om.get_soln('lin', 'mu_l', 'Ax'), 14, [t 'lin.mu_l.Ax']);
    t_is(om.soln.lin.mu_u.Ax, om.get_soln('lin', 'mu_u', 'Ax'), 14, [t 'lin.mu_u.Ax']);
else
    t_skip(29, 'no MILP/MIQP solver installed');
end

if have_feature('quadprog') && have_feature('quadprog', 'vnum') == 7.005
    warning(s1.state, diff_alg_warn_id);
end

t_end;

function TorF = have_miqp_solver()
TorF = have_feature('cplex') || have_feature('glpk') || ...
    have_feature('gurobi') || have_feature('intlinprog') || ...
    have_feature('mosek');

