function mpc = t_case_int
% t_case_int - Case data in internal format used to test ext2int and int2ext.

%   MATPOWER

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	0	0	0	1	1	0	345	1	1.1	0.9;
	2	2	0	0	0	0	1	1	0	345	1	1.1	0.9;
	3	2	0	0	0	0	1	1	0	345	1	1.1	0.9;
	4	1	0	0	0	0	1	1	0	345	1	1.1	0.9;
	5	1	90	30	0	0.2	1	1	0	345	1	1.1	0.9;
	6	1	0	0	0	0	1	1	0	345	1	1.1	0.9;
	7	1	100	35	0	0	1	1	0	345	1	1.1	0.9;
	8	1	0	0	0.1	0	1	1	0	345	1	1.1	0.9;
	9	1	125	50	0	0	1	1	0	345	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	3	85	0	300	-300	1	100	1	270	10	0	0	0	0	0	0	0	0	0	0	0;
	2	163	0	300	-300	1	100	1	300	10	0	0	0	0	0	0	0	0	0	0	0;
	1	0	0	300	-300	1	100	1	250	90	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	4	0	0.0576	0	0	250	250	0	0	1	-360	360;
	4	5	0.017	0.092	0.158	0	250	250	0	0	1	-360	360;
	5	6	0.039	0.17	0.358	150	150	150	0	0	1	-360	360;
	3	6	0	0.0586	0	0	300	300	0	0	1	-360	360;
	6	7	0.0119	0.1008	0.209	40	150	150	0	0	1	-360	360;
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
	2	0	0	2	15	0	0	0	0	0	0	0;
	1	0	0	4	0	0	100	2500	200	5500	250	7250;
	1	0	0	4	0	0	100	2000	200	4403.5	270	6363.5;
];

%% bus names
mpc.bus_name = {
            'BUS 1';
            'BUS 2';
            'BUS 30';
            'BUS 4';
            'BUS 5';
            'BUS 6';
            'BUS 7';
            'BUS 8';
            'BUS 9';
};

%% generator unit type (see GENTYPES)
mpc.gentype = {
	'ST';
	'GT';
	'HY';
};

%% generator fuel type (see GENFUELS)
mpc.genfuel = {
	'coal';
	'ng';
	'hydro';
};

mpc.A = [
	1	2	3	4	5	7	8	9	10	11	12	13	14	15	17	18	19	20	21	22	24	25	26	28	29	30;
	2	4	6	8	10	14	16	18	20	22	24	26	28	30	34	36	38	40	42	44	48	50	52	56	58	60;
];

mpc.N = [
	30	29	28	27	26	24	23	22	21	20	19	18	17	16	14	13	12	11	10	9	7	6	5	3	2	1;
	60	58	56	54	52	48	46	44	42	40	38	36	34	32	28	26	24	22	20	18	14	12	10	6	4	2;
];
