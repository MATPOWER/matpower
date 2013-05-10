function t_qps_matpower(quiet)
%T_QPS_MATPOWER  Tests of QPS_MATPOWER QP solvers.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2010-2011 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

if nargin < 1
    quiet = 0;
end

algs = [100 200 250 400 300 500 600 700];
names = {'BPMPD_MEX', 'MIPS', 'sc-MIPS', 'IPOPT', 'linprog/quadprog', 'CPLEX', 'MOSEK', 'Gurobi'};
check = {'bpmpd', [], [], 'ipopt', 'quadprog', 'cplex', 'mosek', 'gurobi'};

n = 36;
t_begin(n*length(algs), quiet);

for k = 1:length(algs)
    if ~isempty(check{k}) && ~have_fcn(check{k})
        t_skip(n, sprintf('%s not installed', names{k}));
    else
        opt = struct('verbose', 0, 'alg', algs(k));
        if strcmp(names{k}, 'MIPS') || strcmp(names{k}, 'sc-MIPS')
            opt.mips_opt.comptol = 1e-8;
        end
%         if strcmp(names{k}, 'quadprog')
%         end
        if strcmp(names{k}, 'CPLEX')
%           alg = 0;        %% default uses barrier method with NaN bug in lower lim multipliers
            alg = 2;        %% use dual simplex
            mpopt = mpoption('CPLEX_LPMETHOD', alg, 'CPLEX_QPMETHOD', min([4 alg]));
            opt.cplex_opt = cplex_options([], mpopt);
        end
        if strcmp(names{k}, 'MOSEK')
%             alg = 5;        %% use dual simplex
            mpopt = mpoption;
%             mpopt = mpoption(mpopt, 'MOSEK_LP_ALG', alg );
            mpopt = mpoption(mpopt, 'MOSEK_GAP_TOL', 1e-10);
            opt.mosek_opt = mosek_options([], mpopt);
        end

        t = sprintf('%s - 3-d LP : ', names{k});
        %% example from 'doc linprog'
        c = [-5; -4; -6];
        A = [1 -1  1;
             3  2  4;
             3  2  0];
        l = [];
        u = [20; 42; 30];
        xmin = [0; 0; 0];
        x0 = [];
        [x, f, s, out, lam] = qps_matpower([], c, A, l, u, xmin, [], [], opt);
        t_is(s, 1, 12, [t 'success']);
        t_is(x, [0; 15; 3], 6, [t 'x']);
        t_is(f, -78, 6, [t 'f']);
        t_is(lam.mu_l, [0;0;0], 9, [t 'lam.mu_l']);
        t_is(lam.mu_u, [0;1.5;0.5], 9, [t 'lam.mu_u']);
        t_is(lam.lower, [1;0;0], 9, [t 'lam.lower']);
        t_is(lam.upper, zeros(size(x)), 9, [t 'lam.upper']);

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
        t_is(lam.mu_u, [28;4;0]/9, 7, [t 'lam.mu_u']);
        t_is(lam.lower, zeros(size(x)), 8, [t 'lam.lower']);
        t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

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
        t_is(lam.lower, [2.24;0;0;1.7667], 4, [t 'lam.lower']);
        t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

        t = sprintf('%s - (struct) constrained 4-d QP : ', names{k});
        p = struct('H', H, 'A', A, 'l', l, 'u', u, 'xmin', xmin, 'x0', x0, 'opt', opt);
        [x, f, s, out, lam] = qps_matpower(p);
        t_is(s, 1, 12, [t 'success']);
        t_is(x, [0; 2.8; 0.2; 0]/3, 5, [t 'x']);
        t_is(f, 3.29/3, 6, [t 'f']);
        t_is(lam.mu_l, [6.58;0]/3, 6, [t 'lam.mu_l']);
        t_is(lam.mu_u, [0;0], 13, [t 'lam.mu_u']);
        t_is(lam.lower, [2.24;0;0;1.7667], 4, [t 'lam.lower']);
        t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

        t = sprintf('%s - infeasible LP : ', names{k});
        p = struct('A', sparse([1 1]), 'c', [1;1], 'u', -1, 'xmin', [0;0], 'opt', opt);
        [x, f, s, out, lam] = qps_matpower(p);
        t_ok(s <= 0, [t 'no success']);
    end
end

t_end;
