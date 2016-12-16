function mpc = ex_case3a
%EX_CASE3A  Three bus example system.
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
	1	125	0	25	-25	1	100	1	200	0	0	0	0	0	0	0	0	250	250	0	0;
	1	125	0	25	-25	1	100	1	200	0	0	0	0	0	0	0	0	250	250	0	0;
	2	200	0	50	-50	1	100	1	500	0	0	0	0	0	0	0	0	600	600	0	0;
	3	-450	0	0	0	1	100	1	0	-450	0	0	0	0	0	0	0	500	500	0	0;
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
mpc.gencost = [
	2	0	0	3	0.1	0	0;
	2	0	0	3	0.1	0	0;
	2	0	0	3	0.1	0	0;
	2	0	0	3	0	1000	0;
];

%%-----  Reserve Data  -----%%
%% reserve zones, element i, j is 1 if gen j is in zone i, 0 otherwise
mpc.reserves.zones = [
	1	1	1	0;
];

%% reserve requirements for each zone in MW
mpc.reserves.req   = 150;

%% reserve costs in $/MW for each gen that belongs to at least 1 zone
%% (same order as gens, but skipping any gen that does not belong to any zone)
% mpc.reserves.cost  = [	5;	5;	21;	];
% mpc.reserves.cost  = [	5;	5;	16.25;	];
% mpc.reserves.cost  = [	0;	0;	11.25;	];
mpc.reserves.cost  = [	1;	3;	5;	];

%% OPTIONAL max reserve quantities for each gen that belongs to at least 1 zone
%% (same order as gens, but skipping any gen that does not belong to any zone)
mpc.reserves.qty   = [	100;	100;	200;	];
