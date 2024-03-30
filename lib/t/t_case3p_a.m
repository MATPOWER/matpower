function mpc = t_case3p_a
% t_case3p_a - Four bus, unbalanced 3-phase test case.
%
% This data comes from ``4Bus-YY-UnB.DSS``, a modified version (with unbalanced
% load) of ``4Bus-YY-Bal.DSS`` [1], the OpenDSS 4 bus IEEE test case with
% grounded-wye to grounded-wye transformer.
%
% [1] https://sourceforge.net/p/electricdss/code/HEAD/tree/trunk/Distrib/IEEETestCases/4Bus-YY-Bal/4Bus-YY-Bal.DSS

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

mpc.bus = [];
mpc.gen = [];
mpc.branch = [];
mpc.gencost = [];


%%-----  3 Phase Model Data  -----%%
%% system data
mpc.freq = 60;      %% frequency, Hz
mpc.basekVA = 1000; %% system kVA base

%% bus data
%	busid	type	basekV	Vm1	Vm2	Vm3	Va1	Va2	Va3
mpc.bus3p = [
	1	3	12.47	1	1	1	0	-120	120;
	2	1	12.47	1	1	1	0	-120	120;
	3	1	4.16	1	1	1	0	-120	120;
	4	1	4.16	1	1	1	0	-120	120;
];

%% branch data
%	brid	fbus	tbus	status	lcid	len
mpc.line3p = [
	1	1	2	1	1	2000/5280;
	2	3	4	1	1	2500/5280;
];

%% transformer
%	xfid	fbus	tbus	status	R	X	basekVA	basekV
mpc.xfmr3p = [
	1	2	3	1	0.01	0.06	6000	12.47;
];

%% load
%	ldid	ldbus	status	Pd1	Pd2	Pd3	ldpf1	ldpf2	ldpf3
mpc.load3p = [
	1	4	1	1275	1800	2375	0.85	0.9	0.95;
];

%% gen
%	genid	gbus	status	Vg1	Vg2	Vg3	Pg1	Pg2	Pg3	Qg1	Qg2	Qg3
mpc.gen3p = [
	1	1	1	1	1	1	2000	2000	2000	0	0	0;
];

%% line construction
%	lcid	R11	R21	R31	R22	R32	R33	X11	X21	X31	X22	X32	X33	C11	C21	C31	C22	C32	C33
mpc.lc = [
	1	0.457541	0.15594 	0.153474	0.466617	0.157996	0.461462	1.078	0.501648	0.384909	1.04813	0.423624	1.06502	15.0671	-4.86241	-1.85323	15.875	-3.09098	14.3254
];
