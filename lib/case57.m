function [baseMVA, bus, gen, branch, areas, gencost] = case57
%CASE57    Power flow data for IEEE 57 bus test case.
%   Please see 'help caseformat' for details on the case file format.
%   This data was converted from IEEE Common Data Format
%   (ieee57cdf.txt) on 20-Sep-2004 by cdf2matp, rev. 1.11
%   See end of file for warnings generated during conversion.
%
%   Converted from IEEE CDF file from:
%       http://www.ee.washington.edu/research/pstca/
%
%   Manually modified Qmax, Qmin on generator 1 to 200, -140, respectively.
% 
%  08/25/93 UW ARCHIVE           100.0  1961 W IEEE 57 Bus Test Case

%   MATPOWER
%   $Id$

%%-----  Power Flow Data  -----%%
%% system MVA base
baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
bus = [
	1	3	55	17	0	0	1	1.04	0	0	1	1.06	0.94;
	2	2	3	88	0	0	1	1.01	-1.18	0	1	1.06	0.94;
	3	2	41	21	0	0	1	0.985	-5.97	0	1	1.06	0.94;
	4	1	0	0	0	0	1	0.981	-7.32	0	1	1.06	0.94;
	5	1	13	4	0	0	1	0.976	-8.52	0	1	1.06	0.94;
	6	2	75	2	0	0	1	0.98	-8.65	0	1	1.06	0.94;
	7	1	0	0	0	0	1	0.984	-7.58	0	1	1.06	0.94;
	8	2	150	22	0	0	1	1.005	-4.45	0	1	1.06	0.94;
	9	2	121	26	0	0	1	0.98	-9.56	0	1	1.06	0.94;
	10	1	5	2	0	0	1	0.986	-11.43	0	1	1.06	0.94;
	11	1	0	0	0	0	1	0.974	-10.17	0	1	1.06	0.94;
	12	2	377	24	0	0	1	1.015	-10.46	0	1	1.06	0.94;
	13	1	18	2.3	0	0	1	0.979	-9.79	0	1	1.06	0.94;
	14	1	10.5	5.3	0	0	1	0.97	-9.33	0	1	1.06	0.94;
	15	1	22	5	0	0	1	0.988	-7.18	0	1	1.06	0.94;
	16	1	43	3	0	0	1	1.013	-8.85	0	1	1.06	0.94;
	17	1	42	8	0	0	1	1.017	-5.39	0	1	1.06	0.94;
	18	1	27.2	9.8	0	10	1	1.001	-11.71	0	1	1.06	0.94;
	19	1	3.3	0.6	0	0	1	0.97	-13.2	0	1	1.06	0.94;
	20	1	2.3	1	0	0	1	0.964	-13.41	0	1	1.06	0.94;
	21	1	0	0	0	0	1	1.008	-12.89	0	1	1.06	0.94;
	22	1	0	0	0	0	1	1.01	-12.84	0	1	1.06	0.94;
	23	1	6.3	2.1	0	0	1	1.008	-12.91	0	1	1.06	0.94;
	24	1	0	0	0	0	1	0.999	-13.25	0	1	1.06	0.94;
	25	1	6.3	3.2	0	5.9	1	0.982	-18.13	0	1	1.06	0.94;
	26	1	0	0	0	0	1	0.959	-12.95	0	1	1.06	0.94;
	27	1	9.3	0.5	0	0	1	0.982	-11.48	0	1	1.06	0.94;
	28	1	4.6	2.3	0	0	1	0.997	-10.45	0	1	1.06	0.94;
	29	1	17	2.6	0	0	1	1.01	-9.75	0	1	1.06	0.94;
	30	1	3.6	1.8	0	0	1	0.962	-18.68	0	1	1.06	0.94;
	31	1	5.8	2.9	0	0	1	0.936	-19.34	0	1	1.06	0.94;
	32	1	1.6	0.8	0	0	1	0.949	-18.46	0	1	1.06	0.94;
	33	1	3.8	1.9	0	0	1	0.947	-18.5	0	1	1.06	0.94;
	34	1	0	0	0	0	1	0.959	-14.1	0	1	1.06	0.94;
	35	1	6	3	0	0	1	0.966	-13.86	0	1	1.06	0.94;
	36	1	0	0	0	0	1	0.976	-13.59	0	1	1.06	0.94;
	37	1	0	0	0	0	1	0.985	-13.41	0	1	1.06	0.94;
	38	1	14	7	0	0	1	1.013	-12.71	0	1	1.06	0.94;
	39	1	0	0	0	0	1	0.983	-13.46	0	1	1.06	0.94;
	40	1	0	0	0	0	1	0.973	-13.62	0	1	1.06	0.94;
	41	1	6.3	3	0	0	1	0.996	-14.05	0	1	1.06	0.94;
	42	1	7.1	4.4	0	0	1	0.966	-15.5	0	1	1.06	0.94;
	43	1	2	1	0	0	1	1.01	-11.33	0	1	1.06	0.94;
	44	1	12	1.8	0	0	1	1.017	-11.86	0	1	1.06	0.94;
	45	1	0	0	0	0	1	1.036	-9.25	0	1	1.06	0.94;
	46	1	0	0	0	0	1	1.05	-11.89	0	1	1.06	0.94;
	47	1	29.7	11.6	0	0	1	1.033	-12.49	0	1	1.06	0.94;
	48	1	0	0	0	0	1	1.027	-12.59	0	1	1.06	0.94;
	49	1	18	8.5	0	0	1	1.036	-12.92	0	1	1.06	0.94;
	50	1	21	10.5	0	0	1	1.023	-13.39	0	1	1.06	0.94;
	51	1	18	5.3	0	0	1	1.052	-12.52	0	1	1.06	0.94;
	52	1	4.9	2.2	0	0	1	0.98	-11.47	0	1	1.06	0.94;
	53	1	20	10	0	6.3	1	0.971	-12.23	0	1	1.06	0.94;
	54	1	4.1	1.4	0	0	1	0.996	-11.69	0	1	1.06	0.94;
	55	1	6.8	3.4	0	0	1	1.031	-10.78	0	1	1.06	0.94;
	56	1	7.6	2.2	0	0	1	0.968	-16.04	0	1	1.06	0.94;
	57	1	6.7	2	0	0	1	0.965	-16.56	0	1	1.06	0.94;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin
gen = [
	1	128.9	-16.1	200	-140	1.04	100	1	575.88	0;
	2	0	-0.8	50	-17	1.01	100	1	100	0;
	3	40	-1	60	-10	0.985	100	1	140	0;
	6	0	0.8	25	-8	0.98	100	1	100	0;
	8	450	62.1	200	-140	1.005	100	1	550	0;
	9	0	2.2	9	-3	0.98	100	1	100	0;
	12	310	128.5	155	-150	1.015	100	1	410	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status
branch = [
	1	2	0.0083	0.028	0.129	9900	0	0	0	0	1;
	2	3	0.0298	0.085	0.0818	9900	0	0	0	0	1;
	3	4	0.0112	0.0366	0.038	9900	0	0	0	0	1;
	4	5	0.0625	0.132	0.0258	9900	0	0	0	0	1;
	4	6	0.043	0.148	0.0348	9900	0	0	0	0	1;
	6	7	0.02	0.102	0.0276	9900	0	0	0	0	1;
	6	8	0.0339	0.173	0.047	9900	0	0	0	0	1;
	8	9	0.0099	0.0505	0.0548	9900	0	0	0	0	1;
	9	10	0.0369	0.1679	0.044	9900	0	0	0	0	1;
	9	11	0.0258	0.0848	0.0218	9900	0	0	0	0	1;
	9	12	0.0648	0.295	0.0772	9900	0	0	0	0	1;
	9	13	0.0481	0.158	0.0406	9900	0	0	0	0	1;
	13	14	0.0132	0.0434	0.011	9900	0	0	0	0	1;
	13	15	0.0269	0.0869	0.023	9900	0	0	0	0	1;
	1	15	0.0178	0.091	0.0988	9900	0	0	0	0	1;
	1	16	0.0454	0.206	0.0546	9900	0	0	0	0	1;
	1	17	0.0238	0.108	0.0286	9900	0	0	0	0	1;
	3	15	0.0162	0.053	0.0544	9900	0	0	0	0	1;
	4	18	0	0.555	0	9900	0	0	0.97	0	1;
	4	18	0	0.43	0	9900	0	0	0.978	0	1;
	5	6	0.0302	0.0641	0.0124	9900	0	0	0	0	1;
	7	8	0.0139	0.0712	0.0194	9900	0	0	0	0	1;
	10	12	0.0277	0.1262	0.0328	9900	0	0	0	0	1;
	11	13	0.0223	0.0732	0.0188	9900	0	0	0	0	1;
	12	13	0.0178	0.058	0.0604	9900	0	0	0	0	1;
	12	16	0.018	0.0813	0.0216	9900	0	0	0	0	1;
	12	17	0.0397	0.179	0.0476	9900	0	0	0	0	1;
	14	15	0.0171	0.0547	0.0148	9900	0	0	0	0	1;
	18	19	0.461	0.685	0	9900	0	0	0	0	1;
	19	20	0.283	0.434	0	9900	0	0	0	0	1;
	21	20	0	0.7767	0	9900	0	0	1.043	0	1;
	21	22	0.0736	0.117	0	9900	0	0	0	0	1;
	22	23	0.0099	0.0152	0	9900	0	0	0	0	1;
	23	24	0.166	0.256	0.0084	9900	0	0	0	0	1;
	24	25	0	1.182	0	9900	0	0	1	0	1;
	24	25	0	1.23	0	9900	0	0	1	0	1;
	24	26	0	0.0473	0	9900	0	0	1.043	0	1;
	26	27	0.165	0.254	0	9900	0	0	0	0	1;
	27	28	0.0618	0.0954	0	9900	0	0	0	0	1;
	28	29	0.0418	0.0587	0	9900	0	0	0	0	1;
	7	29	0	0.0648	0	9900	0	0	0.967	0	1;
	25	30	0.135	0.202	0	9900	0	0	0	0	1;
	30	31	0.326	0.497	0	9900	0	0	0	0	1;
	31	32	0.507	0.755	0	9900	0	0	0	0	1;
	32	33	0.0392	0.036	0	9900	0	0	0	0	1;
	34	32	0	0.953	0	9900	0	0	0.975	0	1;
	34	35	0.052	0.078	0.0032	9900	0	0	0	0	1;
	35	36	0.043	0.0537	0.0016	9900	0	0	0	0	1;
	36	37	0.029	0.0366	0	9900	0	0	0	0	1;
	37	38	0.0651	0.1009	0.002	9900	0	0	0	0	1;
	37	39	0.0239	0.0379	0	9900	0	0	0	0	1;
	36	40	0.03	0.0466	0	9900	0	0	0	0	1;
	22	38	0.0192	0.0295	0	9900	0	0	0	0	1;
	11	41	0	0.749	0	9900	0	0	0.955	0	1;
	41	42	0.207	0.352	0	9900	0	0	0	0	1;
	41	43	0	0.412	0	9900	0	0	0	0	1;
	38	44	0.0289	0.0585	0.002	9900	0	0	0	0	1;
	15	45	0	0.1042	0	9900	0	0	0.955	0	1;
	14	46	0	0.0735	0	9900	0	0	0.9	0	1;
	46	47	0.023	0.068	0.0032	9900	0	0	0	0	1;
	47	48	0.0182	0.0233	0	9900	0	0	0	0	1;
	48	49	0.0834	0.129	0.0048	9900	0	0	0	0	1;
	49	50	0.0801	0.128	0	9900	0	0	0	0	1;
	50	51	0.1386	0.22	0	9900	0	0	0	0	1;
	10	51	0	0.0712	0	9900	0	0	0.93	0	1;
	13	49	0	0.191	0	9900	0	0	0.895	0	1;
	29	52	0.1442	0.187	0	9900	0	0	0	0	1;
	52	53	0.0762	0.0984	0	9900	0	0	0	0	1;
	53	54	0.1878	0.232	0	9900	0	0	0	0	1;
	54	55	0.1732	0.2265	0	9900	0	0	0	0	1;
	11	43	0	0.153	0	9900	0	0	0.958	0	1;
	44	45	0.0624	0.1242	0.004	9900	0	0	0	0	1;
	40	56	0	1.195	0	9900	0	0	0.958	0	1;
	56	41	0.553	0.549	0	9900	0	0	0	0	1;
	56	42	0.2125	0.354	0	9900	0	0	0	0	1;
	39	57	0	1.355	0	9900	0	0	0.98	0	1;
	57	56	0.174	0.26	0	9900	0	0	0	0	1;
	38	49	0.115	0.177	0.003	9900	0	0	0	0	1;
	38	48	0.0312	0.0482	0	9900	0	0	0	0	1;
	9	55	0	0.1205	0	9900	0	0	0.94	0	1;
];

%%-----  OPF Data  -----%%
%% area data
areas = [
	1	1;
];

%% generator cost data
%	1	startup	shutdown	n	x0	y0	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
gencost = [
	2	0	0	3	0.0775795	20	0;
	2	0	0	3	0.01	40	0;
	2	0	0	3	0.25	20	0;
	2	0	0	3	0.01	40	0;
	2	0	0	3	0.0222222	20	0;
	2	0	0	3	0.01	40	0;
	2	0	0	3	0.0322581	20	0;
];

return;

% Warnings from cdf2matp conversion:
%
% ***** Qmax = Qmin at generator at bus    1 (Qmax set to Qmin + 10)
% ***** area data conversion not yet implemented (creating dummy area data)
% ***** Insufficient generation, setting Pmax at slack bus (bus 1) to 575.88
% ***** MVA limit of branch 1 - 2 not given, set to 9900
% ***** MVA limit of branch 2 - 3 not given, set to 9900
% ***** MVA limit of branch 3 - 4 not given, set to 9900
% ***** MVA limit of branch 4 - 5 not given, set to 9900
% ***** MVA limit of branch 4 - 6 not given, set to 9900
% ***** MVA limit of branch 6 - 7 not given, set to 9900
% ***** MVA limit of branch 6 - 8 not given, set to 9900
% ***** MVA limit of branch 8 - 9 not given, set to 9900
% ***** MVA limit of branch 9 - 10 not given, set to 9900
% ***** MVA limit of branch 9 - 11 not given, set to 9900
% ***** MVA limit of branch 9 - 12 not given, set to 9900
% ***** MVA limit of branch 9 - 13 not given, set to 9900
% ***** MVA limit of branch 13 - 14 not given, set to 9900
% ***** MVA limit of branch 13 - 15 not given, set to 9900
% ***** MVA limit of branch 1 - 15 not given, set to 9900
% ***** MVA limit of branch 1 - 16 not given, set to 9900
% ***** MVA limit of branch 1 - 17 not given, set to 9900
% ***** MVA limit of branch 3 - 15 not given, set to 9900
% ***** MVA limit of branch 4 - 18 not given, set to 9900
% ***** MVA limit of branch 4 - 18 not given, set to 9900
% ***** MVA limit of branch 5 - 6 not given, set to 9900
% ***** MVA limit of branch 7 - 8 not given, set to 9900
% ***** MVA limit of branch 10 - 12 not given, set to 9900
% ***** MVA limit of branch 11 - 13 not given, set to 9900
% ***** MVA limit of branch 12 - 13 not given, set to 9900
% ***** MVA limit of branch 12 - 16 not given, set to 9900
% ***** MVA limit of branch 12 - 17 not given, set to 9900
% ***** MVA limit of branch 14 - 15 not given, set to 9900
% ***** MVA limit of branch 18 - 19 not given, set to 9900
% ***** MVA limit of branch 19 - 20 not given, set to 9900
% ***** MVA limit of branch 21 - 20 not given, set to 9900
% ***** MVA limit of branch 21 - 22 not given, set to 9900
% ***** MVA limit of branch 22 - 23 not given, set to 9900
% ***** MVA limit of branch 23 - 24 not given, set to 9900
% ***** MVA limit of branch 24 - 25 not given, set to 9900
% ***** MVA limit of branch 24 - 25 not given, set to 9900
% ***** MVA limit of branch 24 - 26 not given, set to 9900
% ***** MVA limit of branch 26 - 27 not given, set to 9900
% ***** MVA limit of branch 27 - 28 not given, set to 9900
% ***** MVA limit of branch 28 - 29 not given, set to 9900
% ***** MVA limit of branch 7 - 29 not given, set to 9900
% ***** MVA limit of branch 25 - 30 not given, set to 9900
% ***** MVA limit of branch 30 - 31 not given, set to 9900
% ***** MVA limit of branch 31 - 32 not given, set to 9900
% ***** MVA limit of branch 32 - 33 not given, set to 9900
% ***** MVA limit of branch 34 - 32 not given, set to 9900
% ***** MVA limit of branch 34 - 35 not given, set to 9900
% ***** MVA limit of branch 35 - 36 not given, set to 9900
% ***** MVA limit of branch 36 - 37 not given, set to 9900
% ***** MVA limit of branch 37 - 38 not given, set to 9900
% ***** MVA limit of branch 37 - 39 not given, set to 9900
% ***** MVA limit of branch 36 - 40 not given, set to 9900
% ***** MVA limit of branch 22 - 38 not given, set to 9900
% ***** MVA limit of branch 11 - 41 not given, set to 9900
% ***** MVA limit of branch 41 - 42 not given, set to 9900
% ***** MVA limit of branch 41 - 43 not given, set to 9900
% ***** MVA limit of branch 38 - 44 not given, set to 9900
% ***** MVA limit of branch 15 - 45 not given, set to 9900
% ***** MVA limit of branch 14 - 46 not given, set to 9900
% ***** MVA limit of branch 46 - 47 not given, set to 9900
% ***** MVA limit of branch 47 - 48 not given, set to 9900
% ***** MVA limit of branch 48 - 49 not given, set to 9900
% ***** MVA limit of branch 49 - 50 not given, set to 9900
% ***** MVA limit of branch 50 - 51 not given, set to 9900
% ***** MVA limit of branch 10 - 51 not given, set to 9900
% ***** MVA limit of branch 13 - 49 not given, set to 9900
% ***** MVA limit of branch 29 - 52 not given, set to 9900
% ***** MVA limit of branch 52 - 53 not given, set to 9900
% ***** MVA limit of branch 53 - 54 not given, set to 9900
% ***** MVA limit of branch 54 - 55 not given, set to 9900
% ***** MVA limit of branch 11 - 43 not given, set to 9900
% ***** MVA limit of branch 44 - 45 not given, set to 9900
% ***** MVA limit of branch 40 - 56 not given, set to 9900
% ***** MVA limit of branch 56 - 41 not given, set to 9900
% ***** MVA limit of branch 56 - 42 not given, set to 9900
% ***** MVA limit of branch 39 - 57 not given, set to 9900
% ***** MVA limit of branch 57 - 56 not given, set to 9900
% ***** MVA limit of branch 38 - 49 not given, set to 9900
% ***** MVA limit of branch 38 - 48 not given, set to 9900
% ***** MVA limit of branch 9 - 55 not given, set to 9900
