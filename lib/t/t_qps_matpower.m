function t_qps_matpower(quiet)
%T_QPS_MATPOWER  Tests of QPS_MATPOWER QP solvers.

%   MATPOWER
%   Copyright (c) 2010-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

algs = {'BPMPD', 'MIPS', 250, 'IPOPT', 'OT', 'CPLEX', 'MOSEK', 'GUROBI', 'CLP', 'GLPK'};
names = {'BPMPD_MEX', 'MIPS', 'sc-MIPS', 'IPOPT', 'linprog/quadprog', 'CPLEX', 'MOSEK', 'Gurobi', 'CLP', 'glpk'};
check = {'bpmpd', [], [], 'ipopt', 'quadprog', 'cplex', 'mosek', 'gurobi', 'clp', 'glpk'};
does_qp = [1 1 1 1 1 1 1 1 1 0];

n = 36;
nqp = 28;
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
        if strcmp(names{k}, 'MIPS') || strcmp(names{k}, 'sc-MIPS')
            opt.mips_opt.comptol = 1e-8;
        end
%         if strcmp(names{k}, 'linprog/quadprog')
%             opt.verbose = 2;
%             opt.linprog_opt.Algorithm = 'interior-point';
%             opt.linprog_opt.Algorithm = 'active-set';
%             opt.linprog_opt.Algorithm = 'simplex';
%             opt.linprog_opt.Algorithm = 'dual-simplex';
%         end
        if strcmp(names{k}, 'CPLEX')
%           alg = 0;        %% default uses barrier method with NaN bug in lower lim multipliers
            alg = 2;        %% use dual simplex
            mpopt = mpoption('cplex.lpmethod', alg, 'cplex.qpmethod', min([4 alg]));
            opt.cplex_opt = cplex_options([], mpopt);
        end
        if strcmp(names{k}, 'MOSEK')
            mpopt = mpoption;
%             sc = mosek_symbcon;
%             alg = sc.MSK_OPTIMIZER_DUAL_SIMPLEX;    %% use dual simplex
%             alg = sc.MSK_OPTIMIZER_INTPNT;          %% use interior point
%             mpopt = mpoption(mpopt, 'mosek.lp_alg', alg );
            mpopt = mpoption(mpopt, 'mosek.gap_tol', 1e-10);
%             mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_INTPNT_TOL_PFEAS', 1e-10);
%             mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_INTPNT_TOL_DFEAS', 1e-10);
%             mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_INTPNT_TOL_INFEAS', 1e-10);
%             mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-10);
            vnum = have_fcn('mosek', 'vnum');
            if vnum >= 8
%                 mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_INTPNT_QO_TOL_PFEAS', 1e-10);
%                 mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_INTPNT_QO_TOL_DFEAS', 1e-10);
%                 mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_INTPNT_QO_TOL_INFEAS', 1e-10);
%                 mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_INTPNT_QO_TOL_MU_RED', 1e-10);
                mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_INTPNT_QO_TOL_REL_GAP', 1e-10);
            end
%             opt.verbose = 3;
            opt.mosek_opt = mosek_options([], mpopt);
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
        [x, f, s, out, lam] = qps_matpower([], c, A, l, u, xmin, [], [], opt);
        t_is(s, 1, 12, [t 'success']);
        t_is(x, [0; 15; 3], 6, [t 'x']);
        t_is(f, -78, 6, [t 'f']);
        t_is(lam.mu_l, [0;1.5;0], 9, [t 'lam.mu_l']);
        t_is(lam.mu_u, [0;0;0.5], 9, [t 'lam.mu_u']);
        if strcmp(algs{k}, 'CLP') && ~have_fcn('opti_clp')
            t_skip(2, [t 'lam.lower/upper : MEXCLP does not return multipliers on var bounds']);
        else
            t_is(lam.lower, [1;0;0], 9, [t 'lam.lower']);
            t_is(lam.upper, zeros(size(x)), 9, [t 'lam.upper']);
        end

        if does_qp(k)
            t = sprintf('%s - unconstrained 3-d quadratic : ', names{k});
            %% from http://www.akiti.ca/QuadProgEx0Constr.html
            H = [5 -2 -1; -2 4 3; -1 3 5];
            c = [2; -35; -47];
            x0 = [0; 0; 0];
            [x, f, s, out, lam] = qps_matpower(H, c, [], [], [], [], [], [], opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [3; 5; 7], 8, [t 'x']);
            t_is(f, -249, 13, [t 'f']);
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
            [x, f, s, out, lam] = qps_matpower(H, c, A, l, u, xmin, [], x0, opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [2; 4]/3, 7, [t 'x']);
            t_is(f, -74/9, 6, [t 'f']);
            t_is(lam.mu_l, [0;0;0], 13, [t 'lam.mu_l']);
            t_is(lam.mu_u, [28;4;0]/9, 4, [t 'lam.mu_u']);
            if strcmp(algs{k}, 'CLP') && ~have_fcn('opti_clp')
                t_skip(2, [t 'lam.lower/upper : MEXCLP does not return multipliers on var bounds']);
            else
                t_is(lam.lower, zeros(size(x)), 7, [t 'lam.lower']);
                t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);
            end

            t = sprintf('%s - constrained 4-d QP : ', names{k});
            %% from http://www.jmu.edu/docs/sasdoc/sashtml/iml/chap8/sect12.htm
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
            [x, f, s, out, lam] = qps_matpower(H, c, A, l, u, xmin, [], x0, opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [0; 2.8; 0.2; 0]/3, 5, [t 'x']);
            t_is(f, 3.29/3, 6, [t 'f']);
            t_is(lam.mu_l, [6.58;0]/3, 6, [t 'lam.mu_l']);
            t_is(lam.mu_u, [0;0], 13, [t 'lam.mu_u']);
            if strcmp(algs{k}, 'CLP') && ~have_fcn('opti_clp')
                t_skip(2, [t 'lam.lower/upper : MEXCLP does not return multipliers on var bounds']);
            else
                t_is(lam.lower, [2.24;0;0;1.7667], 4, [t 'lam.lower']);
                t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);
            end

            t = sprintf('%s - (struct) constrained 4-d QP : ', names{k});
            p = struct('H', H, 'A', A, 'l', l, 'u', u, 'xmin', xmin, 'x0', x0, 'opt', opt);
            [x, f, s, out, lam] = qps_matpower(p);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [0; 2.8; 0.2; 0]/3, 5, [t 'x']);
            t_is(f, 3.29/3, 6, [t 'f']);
            t_is(lam.mu_l, [6.58;0]/3, 6, [t 'lam.mu_l']);
            t_is(lam.mu_u, [0;0], 13, [t 'lam.mu_u']);
            if strcmp(algs{k}, 'CLP') && ~have_fcn('opti_clp')
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
        [x, f, s, out, lam] = qps_matpower(p);
        t_ok(s <= 0, [t 'no success']);
    end
end

if have_fcn('quadprog') && have_fcn('quadprog', 'vnum') == 7.005
    warning(s1.state, diff_alg_warn_id);
end

t_end;
