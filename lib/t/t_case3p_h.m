function mpc = t_case3p_h
% t_case3p_h -21 bus hybrid test case, 9 single-phase buses, 12 3-phase buses.
%
% Same as t_case3p_g, except the PV hybrid bus has the PV on the 3-phase side.
% Three buses are hybrid buses, one REF-PQ, one PQ-PV and the other PQ-PQ.
% Solutions of three-phase portions should match t_case3p_a.

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	0	0	0	1	1	0	12.47	1	1.1	0.9;
	2	2	0	0	0	0	1	1	0	12.47	1	1.1	0.9;
	3	1	0	0	0	0	1	1	0	12.47	1	1.1	0.9;
	4	1	0	0	0	0	1	1	0	12.47	1	1.1	0.9;
	5	1	90	30	0	0	1	1	0	12.47	1	1.1	0.9;
	6	1	0	0	0	0	1	1	0	12.47	1	1.1	0.9;
	7	1	100	35	0	0	1	1	0	12.47	1	1.1	0.9;
	8	1	0	0	0	0	1	1	0	12.47	1	1.1	0.9;
	9	1	125	50	0	0	1	1	0	12.47	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	78	27.03	300	-300	1.0	100	1	250	10	0	0	0	0	0	0	0	0	0	0	0;
	2	163	6.54	300	-300	1.007807804819155	100	1	300	10	0	0	0	0	0	0	0	0	0	0	0;
	3	85	-10.95	300	-300	1.0	100	1	270	10	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	4	0	0.0576	0	250	250	250	0	0	1	-360	360;
	4	5	0.017	0.092	0.158	250	250	250	0	0	1	-360	360;
	5	6	0.039	0.17	0.358	150	150	150	0	0	1	-360	360;
	3	6	0	0.0586	0	300	300	300	0	0	1	-360	360;
	6	7	0.0119	0.1008	0.209	150	150	150	0	0	1	-360	360;
	7	8	0.0085	0.072	0.149	250	250	250	0	0	1	-360	360;
	8	2	0	0.0625	0	250	250	250	0	0	1	-360	360;
	8	9	0.032	0.161	0.306	250	250	250	0	0	1	-360	360;
	9	4	0.01	0.085	0.176	250	250	250	0	0	1	-360	360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	1500	0	3	0.11	5	150;
	2	2000	0	3	0.085	1.2	600;
	2	3000	0	3	0.1225	1	335;
];


%%-----  3 Phase Model Data  -----%%
%% system data
mpc.freq = 60;      %% frequency, Hz
mpc.basekVA = 1000; %% system kVA base

%% bus data
%	busid	type	basekV	Vm1	Vm2	Vm3	Va1	Va2	Va3
mpc.bus3p = [
	1	1	12.47	1	1	1	0	-120	120;
	2	1	12.47	1	1	1	0	-120	120;
	3	1	4.16	1	1	1	0	-120	120;
	4	1	4.16	1	1	1	0	-120	120;
	5	2	12.47	1	1	1	0	-120	120;
	6	1	12.47	1	1	1	0	-120	120;
	7	1	4.16	1	1	1	0	-120	120;
	8	1	4.16	1	1	1	0	-120	120;
	9	1	12.47	1	1	1	0	-120	120;
	10	1	12.47	1	1	1	0	-120	120;
	11	1	4.16	1	1	1	0	-120	120;
	12	1	4.16	1	1	1	0	-120	120;
];

%% buslink data
%	linkid	busid	bus3pid	status
mpc.buslink = [
	1	1	1	1;
	2	3	5	1;
	3	8	9	1;
];

%% branch data
%	brid	fbus	tbus	status	lcid	len
mpc.line3p = [
	1	1	2	1	1	2000/5280;
	2	3	4	1	1	2500/5280;
	3	5	6	1	1	2000/5280;
	4	7	8	1	1	2500/5280;
	5	9	10	1	1	2000/5280;
	6	11	12	1	1	2500/5280;
];

%% transformer
%	xfid	fbus	tbus	status	R	X	basekVA	basekV
mpc.xfmr3p = [
	1	2	3	1	0.01	0.06	6000	12.47;
	2	6	7	1	0.01	0.06	6000	12.47;
	3	10	11	1	0.01	0.06	6000	12.47;
];

%% load
%	ldid	ldbus	status	Pd1	Pd2	Pd3	ldpf1	ldpf2	ldpf3
mpc.load3p = [
	1	4	1	1275	1800	2375	0.85	0.9	0.95;
	2	8	1	1275	1800	2375	0.85	0.9	0.95;
	3	12	1	1275	1800	2375	0.85	0.9	0.95;
];

%% gen
%	genid	gbus	status	Vg1	Vg2	Vg3	Pg1	Pg2	Pg3	Qg1	Qg2	Qg3
mpc.gen3p = [
	1	5	1	1	1	1	2000	2000	2000	1	1	1;
];

%% line construction
%	lcid	R11	R21	R31	R22	R32	R33	X11	X21	X31	X22	X32	X33	C11	C21	C31	C22	C32	C33
mpc.lc = [
	1	0.457541	0.15594 	0.153474	0.466617	0.157996	0.461462	1.078	0.501648	0.384909	1.04813	0.423624	1.06502	15.0671	-4.86241	-1.85323	15.875	-3.09098	14.3254
];
