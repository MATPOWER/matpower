function t_modcost(quiet)
%T_MODCOST  Tests for code in MODCOST.

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

n_tests = 162;

t_begin(n_tests, quiet);

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
gencost0 = [
	2	0	0	3	0.01	0.1	1	0	0	0	0	0;
	2	0	0	5	0.0006	0.005	0.04	0.3	2	0	0	0;
	1	0	0	4	0	0	10	200	20	600	30	1200;
	1	0	0	4	-30	-2400	-20	-1800	-10	-1000	0	0;
];

%%-----  scalar values for alpha  -----
gencost = modcost(gencost0, 5, 'SCALE_F');

t = 'modcost SCALE_F - quadratic';
t_is(totcost(gencost, [0;0;0;0])/5, [1;2;0;0], 8, t);
t_is(totcost(gencost, [1;0;0;0])/5, [1.11;2;0;0], 8, t);
t_is(totcost(gencost, [2;0;0;0])/5, [1.24;2;0;0], 8, t);

t = 'modcost SCALE_F - 4th order polynomial';
t_is(totcost(gencost, [0;0;0;0])/5, [1;2;     0;0], 8, t);
t_is(totcost(gencost, [0;1;0;0])/5, [1;2.3456;0;0], 8, t);
t_is(totcost(gencost, [0;2;0;0])/5, [1;2.8096;0;0], 8, t);

t = 'modcost SCALE_F - pwl (gen)';
t_is(totcost(gencost, [0;0;5;0 ])/5, [1;2;100;0], 8, t);
t_is(totcost(gencost, [0;0;10;0])/5, [1;2;200;0], 8, t);
t_is(totcost(gencost, [0;0;15;0])/5, [1;2;400;0], 8, t);
t_is(totcost(gencost, [0;0;20;0])/5, [1;2;600;0], 8, t);
t_is(totcost(gencost, [0;0;25;0])/5, [1;2;900;0], 8, t);
t_is(totcost(gencost, [0;0;30;0])/5, [1;2;1200;0], 8, t);
t_is(totcost(gencost, [0;0;35;0])/5, [1;2;1500;0], 8, t);

t = 'modcost SCALE_F - pwl (load)';
t_is(totcost(gencost, [0;0;0;-5 ])/5, [1;2;0;-500], 8, t);
t_is(totcost(gencost, [0;0;0;-10])/5, [1;2;0;-1000], 8, t);
t_is(totcost(gencost, [0;0;0;-15])/5, [1;2;0;-1400], 8, t);
t_is(totcost(gencost, [0;0;0;-20])/5, [1;2;0;-1800], 8, t);
t_is(totcost(gencost, [0;0;0;-25])/5, [1;2;0;-2100], 8, t);
t_is(totcost(gencost, [0;0;0;-30])/5, [1;2;0;-2400], 8, t);
t_is(totcost(gencost, [0;0;0;-35])/5, [1;2;0;-2700], 8, t);


gencost = modcost(gencost0, 2, 'SCALE_X');

t = 'modcost SCALE_X - quadratic';
t_is(totcost(gencost, [0;0;0;0]*2), [1;2;0;0], 8, t);
t_is(totcost(gencost, [1;0;0;0]*2), [1.11;2;0;0], 8, t);
t_is(totcost(gencost, [2;0;0;0]*2), [1.24;2;0;0], 8, t);

t = 'modcost SCALE_X - 4th order polynomial';
t_is(totcost(gencost, [0;0;0;0]*2), [1;2;     0;0], 8, t);
t_is(totcost(gencost, [0;1;0;0]*2), [1;2.3456;0;0], 8, t);
t_is(totcost(gencost, [0;2;0;0]*2), [1;2.8096;0;0], 8, t);

t = 'modcost SCALE_X - pwl (gen)';
t_is(totcost(gencost, [0;0;5;0 ]*2), [1;2;100;0], 8, t);
t_is(totcost(gencost, [0;0;10;0]*2), [1;2;200;0], 8, t);
t_is(totcost(gencost, [0;0;15;0]*2), [1;2;400;0], 8, t);
t_is(totcost(gencost, [0;0;20;0]*2), [1;2;600;0], 8, t);
t_is(totcost(gencost, [0;0;25;0]*2), [1;2;900;0], 8, t);
t_is(totcost(gencost, [0;0;30;0]*2), [1;2;1200;0], 8, t);
t_is(totcost(gencost, [0;0;35;0]*2), [1;2;1500;0], 8, t);

t = 'modcost SCALE_X - pwl (load)';
t_is(totcost(gencost, [0;0;0;-5 ]*2), [1;2;0;-500], 8, t);
t_is(totcost(gencost, [0;0;0;-10]*2), [1;2;0;-1000], 8, t);
t_is(totcost(gencost, [0;0;0;-15]*2), [1;2;0;-1400], 8, t);
t_is(totcost(gencost, [0;0;0;-20]*2), [1;2;0;-1800], 8, t);
t_is(totcost(gencost, [0;0;0;-25]*2), [1;2;0;-2100], 8, t);
t_is(totcost(gencost, [0;0;0;-30]*2), [1;2;0;-2400], 8, t);
t_is(totcost(gencost, [0;0;0;-35]*2), [1;2;0;-2700], 8, t);


gencost = modcost(gencost0, 3, 'SHIFT_F');

t = 'modcost SHIFT_F - quadratic';
t_is(totcost(gencost, [0;0;0;0])-3, [1;2;0;0], 8, t);
t_is(totcost(gencost, [1;0;0;0])-3, [1.11;2;0;0], 8, t);
t_is(totcost(gencost, [2;0;0;0])-3, [1.24;2;0;0], 8, t);

t = 'modcost SHIFT_F - 4th order polynomial';
t_is(totcost(gencost, [0;0;0;0])-3, [1;2;     0;0], 8, t);
t_is(totcost(gencost, [0;1;0;0])-3, [1;2.3456;0;0], 8, t);
t_is(totcost(gencost, [0;2;0;0])-3, [1;2.8096;0;0], 8, t);

t = 'modcost SHIFT_F - pwl (gen)';
t_is(totcost(gencost, [0;0;5;0 ])-3, [1;2;100;0], 8, t);
t_is(totcost(gencost, [0;0;10;0])-3, [1;2;200;0], 8, t);
t_is(totcost(gencost, [0;0;15;0])-3, [1;2;400;0], 8, t);
t_is(totcost(gencost, [0;0;20;0])-3, [1;2;600;0], 8, t);
t_is(totcost(gencost, [0;0;25;0])-3, [1;2;900;0], 8, t);
t_is(totcost(gencost, [0;0;30;0])-3, [1;2;1200;0], 8, t);
t_is(totcost(gencost, [0;0;35;0])-3, [1;2;1500;0], 8, t);

t = 'modcost SHIFT_F - pwl (load)';
t_is(totcost(gencost, [0;0;0;-5 ])-3, [1;2;0;-500], 8, t);
t_is(totcost(gencost, [0;0;0;-10])-3, [1;2;0;-1000], 8, t);
t_is(totcost(gencost, [0;0;0;-15])-3, [1;2;0;-1400], 8, t);
t_is(totcost(gencost, [0;0;0;-20])-3, [1;2;0;-1800], 8, t);
t_is(totcost(gencost, [0;0;0;-25])-3, [1;2;0;-2100], 8, t);
t_is(totcost(gencost, [0;0;0;-30])-3, [1;2;0;-2400], 8, t);
t_is(totcost(gencost, [0;0;0;-35])-3, [1;2;0;-2700], 8, t);


gencost = modcost(gencost0, -4, 'SHIFT_X');

t = 'modcost SHIFT_X - quadratic';
t_is(totcost(gencost, [0;0;0;0]-4), [1;2;0;0], 8, t);
t_is(totcost(gencost, [1;0;0;0]-4), [1.11;2;0;0], 8, t);
t_is(totcost(gencost, [2;0;0;0]-4), [1.24;2;0;0], 8, t);

t = 'modcost SHIFT_X - 4th order polynomial';
t_is(totcost(gencost, [0;0;0;0]-4), [1;2;     0;0], 8, t);
t_is(totcost(gencost, [0;1;0;0]-4), [1;2.3456;0;0], 8, t);
t_is(totcost(gencost, [0;2;0;0]-4), [1;2.8096;0;0], 8, t);

t = 'modcost SHIFT_X - pwl (gen)';
t_is(totcost(gencost, [0;0;5;0 ]-4), [1;2;100;0], 8, t);
t_is(totcost(gencost, [0;0;10;0]-4), [1;2;200;0], 8, t);
t_is(totcost(gencost, [0;0;15;0]-4), [1;2;400;0], 8, t);
t_is(totcost(gencost, [0;0;20;0]-4), [1;2;600;0], 8, t);
t_is(totcost(gencost, [0;0;25;0]-4), [1;2;900;0], 8, t);
t_is(totcost(gencost, [0;0;30;0]-4), [1;2;1200;0], 8, t);
t_is(totcost(gencost, [0;0;35;0]-4), [1;2;1500;0], 8, t);

t = 'modcost SHIFT_X - pwl (load)';
t_is(totcost(gencost, [0;0;0;-5 ]-4), [1;2;0;-500], 8, t);
t_is(totcost(gencost, [0;0;0;-10]-4), [1;2;0;-1000], 8, t);
t_is(totcost(gencost, [0;0;0;-15]-4), [1;2;0;-1400], 8, t);
t_is(totcost(gencost, [0;0;0;-20]-4), [1;2;0;-1800], 8, t);
t_is(totcost(gencost, [0;0;0;-25]-4), [1;2;0;-2100], 8, t);
t_is(totcost(gencost, [0;0;0;-30]-4), [1;2;0;-2400], 8, t);
t_is(totcost(gencost, [0;0;0;-35]-4), [1;2;0;-2700], 8, t);

t = 'modcost empty gencost';
gencost = modcost([], 7);
t_ok(isempty(gencost), t);

%%-----  vector values for alpha  -----
alpha = [10; 9; 8; 7];
gencost = modcost(gencost0, alpha, 'SCALE_F');

t = 'modcost vector SCALE_F - quadratic';
t_is(totcost(gencost, [0;0;0;0])./alpha, [1;2;0;0], 8, t);
t_is(totcost(gencost, [1;0;0;0])./alpha, [1.11;2;0;0], 8, t);
t_is(totcost(gencost, [2;0;0;0])./alpha, [1.24;2;0;0], 8, t);

t = 'modcost vector SCALE_F - 4th order polynomial';
t_is(totcost(gencost, [0;0;0;0])./alpha, [1;2;     0;0], 8, t);
t_is(totcost(gencost, [0;1;0;0])./alpha, [1;2.3456;0;0], 8, t);
t_is(totcost(gencost, [0;2;0;0])./alpha, [1;2.8096;0;0], 8, t);

t = 'modcost vector SCALE_F - pwl (gen)';
t_is(totcost(gencost, [0;0;5;0 ])./alpha, [1;2;100;0], 8, t);
t_is(totcost(gencost, [0;0;10;0])./alpha, [1;2;200;0], 8, t);
t_is(totcost(gencost, [0;0;15;0])./alpha, [1;2;400;0], 8, t);
t_is(totcost(gencost, [0;0;20;0])./alpha, [1;2;600;0], 8, t);
t_is(totcost(gencost, [0;0;25;0])./alpha, [1;2;900;0], 8, t);
t_is(totcost(gencost, [0;0;30;0])./alpha, [1;2;1200;0], 8, t);
t_is(totcost(gencost, [0;0;35;0])./alpha, [1;2;1500;0], 8, t);

t = 'modcost vector SCALE_F - pwl (load)';
t_is(totcost(gencost, [0;0;0;-5 ])./alpha, [1;2;0;-500], 8, t);
t_is(totcost(gencost, [0;0;0;-10])./alpha, [1;2;0;-1000], 8, t);
t_is(totcost(gencost, [0;0;0;-15])./alpha, [1;2;0;-1400], 8, t);
t_is(totcost(gencost, [0;0;0;-20])./alpha, [1;2;0;-1800], 8, t);
t_is(totcost(gencost, [0;0;0;-25])./alpha, [1;2;0;-2100], 8, t);
t_is(totcost(gencost, [0;0;0;-30])./alpha, [1;2;0;-2400], 8, t);
t_is(totcost(gencost, [0;0;0;-35])./alpha, [1;2;0;-2700], 8, t);


gencost = modcost(gencost0, alpha, 'SCALE_X');

t = 'modcost vector SCALE_X - quadratic';
t_is(totcost(gencost, [0;0;0;0].*alpha), [1;2;0;0], 8, t);
t_is(totcost(gencost, [1;0;0;0].*alpha), [1.11;2;0;0], 8, t);
t_is(totcost(gencost, [2;0;0;0].*alpha), [1.24;2;0;0], 8, t);

t = 'modcost vector SCALE_X - 4th order polynomial';
t_is(totcost(gencost, [0;0;0;0].*alpha), [1;2;     0;0], 8, t);
t_is(totcost(gencost, [0;1;0;0].*alpha), [1;2.3456;0;0], 8, t);
t_is(totcost(gencost, [0;2;0;0].*alpha), [1;2.8096;0;0], 8, t);

t = 'modcost vector SCALE_X - pwl (gen)';
t_is(totcost(gencost, [0;0;5;0 ].*alpha), [1;2;100;0], 8, t);
t_is(totcost(gencost, [0;0;10;0].*alpha), [1;2;200;0], 8, t);
t_is(totcost(gencost, [0;0;15;0].*alpha), [1;2;400;0], 8, t);
t_is(totcost(gencost, [0;0;20;0].*alpha), [1;2;600;0], 8, t);
t_is(totcost(gencost, [0;0;25;0].*alpha), [1;2;900;0], 8, t);
t_is(totcost(gencost, [0;0;30;0].*alpha), [1;2;1200;0], 8, t);
t_is(totcost(gencost, [0;0;35;0].*alpha), [1;2;1500;0], 8, t);

t = 'modcost vector SCALE_X - pwl (load)';
t_is(totcost(gencost, [0;0;0;-5 ].*alpha), [1;2;0;-500], 8, t);
t_is(totcost(gencost, [0;0;0;-10].*alpha), [1;2;0;-1000], 8, t);
t_is(totcost(gencost, [0;0;0;-15].*alpha), [1;2;0;-1400], 8, t);
t_is(totcost(gencost, [0;0;0;-20].*alpha), [1;2;0;-1800], 8, t);
t_is(totcost(gencost, [0;0;0;-25].*alpha), [1;2;0;-2100], 8, t);
t_is(totcost(gencost, [0;0;0;-30].*alpha), [1;2;0;-2400], 8, t);
t_is(totcost(gencost, [0;0;0;-35].*alpha), [1;2;0;-2700], 8, t);


gencost = modcost(gencost0, alpha, 'SHIFT_F');

t = 'modcost vector SHIFT_F - quadratic';
t_is(totcost(gencost, [0;0;0;0])-alpha, [1;2;0;0], 8, t);
t_is(totcost(gencost, [1;0;0;0])-alpha, [1.11;2;0;0], 8, t);
t_is(totcost(gencost, [2;0;0;0])-alpha, [1.24;2;0;0], 8, t);

t = 'modcost vector SHIFT_F - 4th order polynomial';
t_is(totcost(gencost, [0;0;0;0])-alpha, [1;2;     0;0], 8, t);
t_is(totcost(gencost, [0;1;0;0])-alpha, [1;2.3456;0;0], 8, t);
t_is(totcost(gencost, [0;2;0;0])-alpha, [1;2.8096;0;0], 8, t);

t = 'modcost vector SHIFT_F - pwl (gen)';
t_is(totcost(gencost, [0;0;5;0 ])-alpha, [1;2;100;0], 8, t);
t_is(totcost(gencost, [0;0;10;0])-alpha, [1;2;200;0], 8, t);
t_is(totcost(gencost, [0;0;15;0])-alpha, [1;2;400;0], 8, t);
t_is(totcost(gencost, [0;0;20;0])-alpha, [1;2;600;0], 8, t);
t_is(totcost(gencost, [0;0;25;0])-alpha, [1;2;900;0], 8, t);
t_is(totcost(gencost, [0;0;30;0])-alpha, [1;2;1200;0], 8, t);
t_is(totcost(gencost, [0;0;35;0])-alpha, [1;2;1500;0], 8, t);

t = 'modcost vector SHIFT_F - pwl (load)';
t_is(totcost(gencost, [0;0;0;-5 ])-alpha, [1;2;0;-500], 8, t);
t_is(totcost(gencost, [0;0;0;-10])-alpha, [1;2;0;-1000], 8, t);
t_is(totcost(gencost, [0;0;0;-15])-alpha, [1;2;0;-1400], 8, t);
t_is(totcost(gencost, [0;0;0;-20])-alpha, [1;2;0;-1800], 8, t);
t_is(totcost(gencost, [0;0;0;-25])-alpha, [1;2;0;-2100], 8, t);
t_is(totcost(gencost, [0;0;0;-30])-alpha, [1;2;0;-2400], 8, t);
t_is(totcost(gencost, [0;0;0;-35])-alpha, [1;2;0;-2700], 8, t);


gencost = modcost(gencost0, -alpha, 'SHIFT_X');

t = 'modcost vector SHIFT_X - quadratic';
t_is(totcost(gencost, [0;0;0;0]-alpha), [1;2;0;0], 8, t);
t_is(totcost(gencost, [1;0;0;0]-alpha), [1.11;2;0;0], 8, t);
t_is(totcost(gencost, [2;0;0;0]-alpha), [1.24;2;0;0], 8, t);

t = 'modcost vector SHIFT_X - 4th order polynomial';
t_is(totcost(gencost, [0;0;0;0]-alpha), [1;2;     0;0], 8, t);
t_is(totcost(gencost, [0;1;0;0]-alpha), [1;2.3456;0;0], 8, t);
t_is(totcost(gencost, [0;2;0;0]-alpha), [1;2.8096;0;0], 8, t);

t = 'modcost vector SHIFT_X - pwl (gen)';
t_is(totcost(gencost, [0;0;5;0 ]-alpha), [1;2;100;0], 8, t);
t_is(totcost(gencost, [0;0;10;0]-alpha), [1;2;200;0], 8, t);
t_is(totcost(gencost, [0;0;15;0]-alpha), [1;2;400;0], 8, t);
t_is(totcost(gencost, [0;0;20;0]-alpha), [1;2;600;0], 8, t);
t_is(totcost(gencost, [0;0;25;0]-alpha), [1;2;900;0], 8, t);
t_is(totcost(gencost, [0;0;30;0]-alpha), [1;2;1200;0], 8, t);
t_is(totcost(gencost, [0;0;35;0]-alpha), [1;2;1500;0], 8, t);

t = 'modcost vector SHIFT_X - pwl (load)';
t_is(totcost(gencost, [0;0;0;-5 ]-alpha), [1;2;0;-500], 8, t);
t_is(totcost(gencost, [0;0;0;-10]-alpha), [1;2;0;-1000], 8, t);
t_is(totcost(gencost, [0;0;0;-15]-alpha), [1;2;0;-1400], 8, t);
t_is(totcost(gencost, [0;0;0;-20]-alpha), [1;2;0;-1800], 8, t);
t_is(totcost(gencost, [0;0;0;-25]-alpha), [1;2;0;-2100], 8, t);
t_is(totcost(gencost, [0;0;0;-30]-alpha), [1;2;0;-2400], 8, t);
t_is(totcost(gencost, [0;0;0;-35]-alpha), [1;2;0;-2700], 8, t);

t = 'modcost vector empty gencost';
gencost = modcost([], alpha);
t_ok(isempty(gencost), t);

t_end;
