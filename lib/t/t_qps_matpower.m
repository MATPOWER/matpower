function t_qps_matpower(quiet)
%T_QPS_MATPOWER  Tests of qps_matpower QP solvers

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2010 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

algs = [100 200 250 300];
names = {'BPMPD_MEX', 'MIPS', 'sc-MIPS', 'quadprog'};
check = {'bpmpd', [], [], 'quadprog'};

n = 35;
t_begin(n*length(algs), quiet);

opt.verbose = 0;
opt.mips_opt.comptol = 1e-8;

for k = 1:length(algs)
	if ~isempty(check{k}) && ~have_fcn(check{k})
		t_skip(n, sprintf('%s not installed', names{k}));
	else
		opt.alg = algs(k);
		if strcmp(names{k}, 'quadprog')
			opt.ot_opt = optimset('LargeScale', 'off');
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
		t_ok(s, [t 'success']);
		t_is(x, [0; 15; 3], 7, [t 'x']);
		t_is(f, -78, 7, [t 'f']);
		t_is(lam.mu_l, [0;0;0], 13, [t 'lam.mu_l']);
		t_is(lam.mu_u, [0;1.5;0.5], 9, [t 'lam.mu_u']);
		t_is(lam.lower, [1;0;0], 9, [t 'lam.lower']);
		t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

		t = sprintf('%s - unconstrained 3-d quadratic : ', names{k});
		%% from http://www.akiti.ca/QuadProgEx0Constr.html
		H = [5 -2 -1; -2 4 3; -1 3 5];
		c = [2; -35; -47];
		x0 = [0; 0; 0];
		[x, f, s, out, lam] = qps_matpower(H, c, [], [], [], [], [], [], opt);
		t_ok(s, [t 'success']);
		t_is(x, [3; 5; 7], 8, [t 'x']);
		t_is(f, -249, 13, [t 'f']);
		t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
		t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
		t_is(lam.lower, zeros(size(x)), 13, [t 'lam.lower']);
		t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);
		
		t = sprintf('%s - constrained 2-d QP : ', names{k});
		%% example from 'doc quadprog'
		H = [   1 	-1;
				-1	2 	];
		c = [-2; -6];
		A = [   1	1;
				-1	2;
				2	1	];
		l = [];
		u = [2; 2; 3];
		xmin = [0; 0];
		x0 = [];
		[x, f, s, out, lam] = qps_matpower(H, c, A, l, u, xmin, [], x0, opt);
		t_ok(s, [t 'success']);
		t_is(x, [2; 4]/3, 7, [t 'x']);
		t_is(f, -74/9, 6, [t 'f']);
		t_is(lam.mu_l, [0;0;0], 13, [t 'lam.mu_l']);
		t_is(lam.mu_u, [28;4;0]/9, 7, [t 'lam.mu_u']);
		t_is(lam.lower, zeros(size(x)), 8, [t 'lam.lower']);
		t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

		t = sprintf('%s - constrained 4-d QP : ', names{k});
		%% from http://www.uc.edu/sashtml/iml/chap8/sect12.htm
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
		t_ok(s, [t 'success']);
		t_is(x, [0; 2.8; 0.2; 0]/3, 6, [t 'x']);
		t_is(f, 3.29/3, 6, [t 'f']);
		t_is(lam.mu_l, [6.58;0]/3, 6, [t 'lam.mu_l']);
		t_is(lam.mu_u, [0;0], 13, [t 'lam.mu_u']);
		t_is(lam.lower, [2.24;0;0;1.7667], 4, [t 'lam.lower']);
		t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

		t = sprintf('%s - (struct) constrained 4-d QP : ', names{k});
		p = struct('H', H, 'A', A, 'l', l, 'u', u, 'xmin', xmin, 'x0', x0, 'opt', opt);
		[x, f, s, out, lam] = qps_matpower(p);
		t_ok(s, [t 'success']);
		t_is(x, [0; 2.8; 0.2; 0]/3, 6, [t 'x']);
		t_is(f, 3.29/3, 6, [t 'f']);
		t_is(lam.mu_l, [6.58;0]/3, 6, [t 'lam.mu_l']);
		t_is(lam.mu_u, [0;0], 13, [t 'lam.mu_u']);
		t_is(lam.lower, [2.24;0;0;1.7667], 4, [t 'lam.lower']);
		t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);
	end
end

t_end;
