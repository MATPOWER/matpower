function t_om_solve_miqps(quiet)
%T_OM_SOLVE_MIQPS  Tests of MIQP solvers via OM.SOLVE().

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

algs = {'DEFAULT', 'CPLEX', 'MOSEK', 'GUROBI', 'GLPK', 'OT'};
names = {'DEFAULT', 'CPLEX', 'MOSEK', 'Gurobi', 'glpk', 'intlin/lin/quadprog'};
check = {[], 'cplex', 'mosek', 'gurobi', 'glpk', 'intlinprog'};
does_qp = [0 1 1 1 0 0];
if have_fcn('gurobi') || have_fcn('cplex') || have_fcn('mosek')
    does_qp(1) = 1;
end

n = 12;
nmiqp = 6;
t_begin(n*length(algs), quiet);

diff_alg_warn_id = 'optim:linprog:WillRunDiffAlg';
if have_fcn('quadprog') && have_fcn('quadprog', 'vnum') == 7.005
    s1 = warning('query', diff_alg_warn_id);
    warning('off', diff_alg_warn_id);
end

for k = 1:length(algs)
    if ~isempty(check{k}) && ~have_fcn(check{k})
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
        if strcmp(names{k}, 'CPLEX')
%           alg = 0;        %% default uses barrier method with NaN bug in lower lim multipliers
            alg = 2;        %% use dual simplex
            mpopt.cplex.lpmethod = alg;
            mpopt.cplex.qpmethod = min([4 alg]);
            opt.cplex_opt = cplex_options([], mpopt);
        end
        if strcmp(names{k}, 'MOSEK')
%             sc = mosek_symbcon;
%             alg = sc.MSK_OPTIMIZER_DUAL_SIMPLEX;    %% use dual simplex
%             alg = sc.MSK_OPTIMIZER_INTPNT;          %% use interior point
%             mpopt.mosek.lp_alg = alg;
            mpopt.mosek.gap_tol = 1e-10;
%             mpopt.mosek.opts.MSK_DPAR_INTPNT_TOL_PFEAS = 1e-10;
%             mpopt.mosek.opts.MSK_DPAR_INTPNT_TOL_DFEAS = 1e-10;
%             mpopt.mosek.opts.MSK_DPAR_INTPNT_TOL_INFEAS = 1e-10;
%             mpopt.mosek.opts.MSK_DPAR_INTPNT_TOL_REL_GAP = 1e-10;
            vnum = have_fcn('mosek', 'vnum');
            if vnum >= 8
%                 mpopt.mosek.opts.MSK_DPAR_INTPNT_QO_TOL_PFEAS = 1e-10;
%                 mpopt.mosek.opts.MSK_DPAR_INTPNT_QO_TOL_DFEAS = 1e-10;
%                 mpopt.mosek.opts.MSK_DPAR_INTPNT_QO_TOL_INFEAS = 1e-10;
%                 mpopt.mosek.opts.MSK_DPAR_INTPNT_QO_TOL_MU_RED = 1e-10;
                mpopt.mosek.opts.MSK_DPAR_INTPNT_QO_TOL_REL_GAP = 1e-10;
            end
%             opt.verbose = 3;
            opt.mosek_opt = mosek_options([], mpopt);
        end

% opt.verbose = 3;
        t = sprintf('%s - 2-d MILP : ', names{k});
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
        t_is(lam.mu_l, [0; 0], 12, [t 'lam.mu_l']);
        t_is(lam.mu_u, [0; 0], 12, [t 'lam.mu_u']);
        t_is(lam.lower, [0; 0], 12, [t 'lam.lower']);
        t_is(lam.upper, [2; 3], 12, [t 'lam.upper']);

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
            t_is(lam.mu_l, [466; 0; 0], 6, [t 'lam.mu_l']);
            t_is(lam.mu_u, [0; 272; 0], 6, [t 'lam.mu_u']);
            t_is(lam.lower, [0; 0; 349.5; 4350], 5, [t 'lam.lower']);
            t_is(lam.upper, [0; 0; 0; 0], 7, [t 'lam.upper']);
        else
            t_skip(nmiqp, sprintf('%s does not handle MIQP problems', names{k}));
        end
% opt.verbose = 0;
    end
end

if have_fcn('quadprog') && have_fcn('quadprog', 'vnum') == 7.005
    warning(s1.state, diff_alg_warn_id);
end

t_end;
