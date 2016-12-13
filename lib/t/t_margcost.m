function t_margcost(quiet)
%T_MARGCOST  Tests for code in MARGCOST.

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

n_tests = 22;

t_begin(n_tests, quiet);

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
gencost = [
	2	0	0	3	0.01	0.1	1	0	0	0	0	0;
	2	0	0	5	0.0006	0.005	0.04	0.3	2	0	0	0;
	1	0	0	4	0	0	10	200	20	600	30	1200;
	1	0	0	4	-30	-2400	-20	-1800	-10	-1000	0	0;
];

t = 'margcost - quadratic';
t_is(margcost(gencost, [0;0;0;0]), [0.1;0.3;20;100], 8, t);
t_is(margcost(gencost, [1;0;0;0]), [0.12;0.3;20;100], 8, t);
t_is(margcost(gencost, [2;0;0;0]), [0.14;0.3;20;100], 8, t);

t = 'margcost - 4th order polynomial';
t_is(margcost(gencost, [0;0;0;0]), [0.1;0.3;20;100], 8, t);
t_is(margcost(gencost, [0;1;0;0]), [0.1;0.3974;20;100], 8, t);
t_is(margcost(gencost, [0;2;0;0]), [0.1;0.5392;20;100], 8, t);

t = 'margcost - pwl (gen)';
t_is(margcost(gencost, [0;0;-10;0 ]), [0.1;0.3;20;100], 8, t);
t_is(margcost(gencost, [0;0;5;0 ]),   [0.1;0.3;20;100], 8, t);
t_is(margcost(gencost, [0;0;10;0]),   [0.1;0.3;20;100], 8, t);
t_is(margcost(gencost, [0;0;15;0]),   [0.1;0.3;40;100], 8, t);
t_is(margcost(gencost, [0;0;20;0]),   [0.1;0.3;40;100], 8, t);
t_is(margcost(gencost, [0;0;25;0]),   [0.1;0.3;60;100], 8, t);
t_is(margcost(gencost, [0;0;30;0]),   [0.1;0.3;60;100], 8, t);
t_is(margcost(gencost, [0;0;35;0]),   [0.1;0.3;60;100], 8, t);

t = 'margcost - pwl (load)';
t_is(margcost(gencost, [0;0;0;10 ]), [0.1;0.3;20;100], 8, t);
t_is(margcost(gencost, [0;0;0;-5 ]), [0.1;0.3;20;100], 8, t);
t_is(margcost(gencost, [0;0;0;-10]), [0.1;0.3;20;80], 8, t);
t_is(margcost(gencost, [0;0;0;-15]), [0.1;0.3;20;80], 8, t);
t_is(margcost(gencost, [0;0;0;-20]), [0.1;0.3;20;60], 8, t);
t_is(margcost(gencost, [0;0;0;-25]), [0.1;0.3;20;60], 8, t);
t_is(margcost(gencost, [0;0;0;-30]), [0.1;0.3;20;60], 8, t);
t_is(margcost(gencost, [0;0;0;-35]), [0.1;0.3;20;60], 8, t);

t_end;
