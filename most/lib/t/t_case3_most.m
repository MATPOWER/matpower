function mpc = t_case3_most
%T_CASE3_MOST  Three bus example system
%   Please see CASEFORMAT for details on the case file format.

%   MOST
%   Copyright (c) 2015-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	0	0	0	1	1	0	135	1	1.05	0.95;
	2	2	0	0	0	0	1	1	0	135	1	1.05	0.95;
	3	2	0	0	0	0	1	1	0	135	1	1.05	0.95;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	250	0	50	-50	1	100	1	400	0	0	0	0	0	0	0	0	100	100	0	0;
	2	200	0	50	-50	1	100	1	300	0	0	0	0	0	0	0	0	100	100	0	0;
	3	-450	0	0	0	1	100	1	0	-450	0	0	0	0	0	0	0	450	450	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.005	0.01	0	300	300	300	0	0	1	-360	360;
	1	3	0.005	0.01	0	240	240	240	0	0	1	-360	360;
	2	3	0.005	0.01	0	300	300	300	0	0	1	-360	360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
% mpc.gencost = [
% 	1	0	0	2	0	0	1	25;
% 	1	0	0	2	0	0	1	30;
% 	1	0	0	2	-1	-10000	0	0;
% ];
mpc.gencost = [
	2	0	0	2	25	0;
	2	0	0	2	30	0;
	2	0	0	2	1000	0;
];
