function t_off2case(quiet)
%T_OFF2CASE  Tests for code in OFF2CASE.

%   MATPOWER
%   Copyright (c) 2005-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

n_tests = 35;

t_begin(n_tests, quiet);

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
gen0 = [
	1	10	0	60	-15	1	100	1	60	10	0	0	0	0	0	0	0	0	0	0	0;
	2	10	0	60	-15	1	100	1	60	12	0	0	0	0	0	0	0	0	0	0	0;
	7	-30	-15	0	-15	1	100	1	0	-30	0	0	0	0	0	0	0	0	0	0	0;
	13	10	0	60	-15	1	100	1	60	12	0	0	0	0	0	0	0	0	0	0	0;
	30	-30	7.5	7.5	0	1	100	1	0	-30	0	0	0	0	0	0	0	0	0	0	0;
];
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
gencost0 = [
	1	0	0	4	0	0	12	240	36	1200	60	2400;
	1	100	0	4	0	0	12	240	36	1200	60	2400;
	1	0	0	4	-30	0	-20	1000	-10	2000	0	3000;
	1	0	0	4	0	0	12	240	36	1200	60	2400;
	1	0	50	4	-30	0	-20	1000	-10	2000	0	3000;
];

if ~have_fcn('smartmarket')
    t_skip(n_tests, 'smartmarket code not available');
else
	t = 'isload()';
	t_is(isload(gen0), [0;0;1;0;1], 8, t);
	
	G = find(~isload(gen0) );
	L = find( isload(gen0) );
	nGL = length(G) + length(L);
	
	t = 'P offers only';
	offers.P.qty = [25; 26; 27];
	offers.P.prc = [10; 50; 100];
	[gen, gencost] = off2case(gen0, gencost0, offers);
	
	gen1 = gen0;
	gen1(G, PMAX) = offers.P.qty;
	gen1(L, GEN_STATUS) = 0;
	t_is( gen, gen1, 8, [t ' - gen'] );
	
	gencost1 = gencost0;
	gencost1(G, NCOST:(NCOST+8)) = [[2 0 0 25 250; 2 0 0 26 1300; 2 0 0 27 2700] zeros(3,4)];
	
	t_is( gencost, gencost1, 8, [t ' - gencost'] );
	
	offers.P.qty = [25; 26; 0; 27; 0];
	offers.P.prc = [10; 50; 0; 100; 0];
	[gen, gencost] = off2case(gen0, gencost0, offers);
	t_is( gen, gen1, 8, [t ' (all rows in offer) - gen'] );
	t_is( gencost, gencost1, 8, [t ' (all rows in offer) - gencost'] );
	
	t = 'P offers only (GEN_STATUS=0 for 0 qty offer)';
	offers.P.qty = [0; 26; 27];
	offers.P.prc = [10; 50; 100];
	[gen, gencost] = off2case(gen0, gencost0, offers);
	
	gen1 = gen0;
	gen1(G(2:3), PMAX) = offers.P.qty(2:3);
	gen1(G(1), GEN_STATUS) = 0;
	gen1(L, GEN_STATUS) = 0;
	t_is( gen, gen1, 8, [t ' - gen'] );
	
	gencost1 = gencost0;
	gencost1(G(2:3), NCOST:(NCOST+8)) = [[2 0 0 26 1300; 2 0 0 27 2700] zeros(2,4)];
	
	t_is( gencost, gencost1, 8, [t ' - gencost'] );

	t = 'P offers, lim.P.max_offer';
	offers.P.qty = [25; 26; 27];
	offers.P.prc = [10; 50; 100];
	lim.P.max_offer = 75;
	[gen, gencost] = off2case(gen0, gencost0, offers, [], lim);
	
	gen1 = gen0;
	gen1(G(1:2), PMAX) = offers.P.qty(1:2, :);
	gen1([G(3); L], GEN_STATUS) = 0;
	t_is( gen, gen1, 8, [t ' - gen'] );
	
	gencost1 = gencost0;
	gencost1(G(1:2), NCOST:(NCOST+8)) = [[2 0 0 25 250; 2 0 0 26 1300] zeros(2,4)];
	t_is( gencost, gencost1, 8, [t ' - gencost'] );
	
	t = 'P offers & P bids';
	bids.P.qty = [20; 28];
	bids.P.prc = [100; 10];
	[gen, gencost] = off2case(gen0, gencost0, offers, bids);
	
	gen1 = gen0;
	gen1(G, PMAX) = offers.P.qty;
	gen1(L, [PMIN QMIN QMAX]) = [-20 -10 0; -28 0 7];
	t_is( gen, gen1, 8, [t ' - gen'] );
	
	gencost1 = gencost0(:, 1:8);
	gencost1(G, NCOST:(NCOST+4)) = [2 0 0 25 250; 2 0 0 26 1300; 2 0 0 27 2700];
	gencost1(L, NCOST:(NCOST+4)) = [2 -20 -2000 0 0; 2 -28 -280 0 0];
	t_is( gencost, gencost1, 8, [t ' - gencost'] );
	
	t = 'P offers & P bids (all rows in bid)';
	bids.P.qty = [0; 0; 20; 0; 28];
	bids.P.prc = [0; 0; 100; 0; 10];
	[gen, gencost] = off2case(gen0, gencost0, offers, bids);
	
	t_is( gen, gen1, 8, [t ' - gen'] );
	t_is( gencost, gencost1, 8, [t ' - gencost'] );
	
	t = 'P offers & P bids (GEN_STATUS=0 for 0 qty bid)';
	bids.P.qty = [0; 28];
	bids.P.prc = [100; 10];
	[gen, gencost] = off2case(gen0, gencost0, offers, bids);
	
	gen1 = gen0;
	gen1(G, PMAX) = offers.P.qty;
	gen1(L(1), GEN_STATUS) = 0;
	gen1(L(2), [PMIN QMIN QMAX]) = [-28 0 7];
	t_is( gen, gen1, 8, [t ' - gen'] );
	
	gencost1 = gencost0;
	gencost1(G, NCOST:(NCOST+8)) = [[2 0 0 25 250; 2 0 0 26 1300; 2 0 0 27 2700] zeros(3,4)];
	gencost1(L(2), NCOST:(NCOST+8)) = [[2 -28 -280 0 0] zeros(1,4)];
	t_is( gencost, gencost1, 8, [t ' - gencost'] );
	
	t = 'P offers & P bids (1 gen with both)';
	gen2 = gen0;
	gen2(2, PMIN) = -5;
	bids.P.qty = [0; 3; 20; 0; 28];
	bids.P.prc = [0; 50; 100; 0; 10];
	[gen, gencost] = off2case(gen2, gencost0, offers, bids);
	
	gen1 = gen2;
	gen1(G, PMAX) = offers.P.qty;
	gen1(2, PMIN) = -sum(bids.P.qty(2, :));
	gen1(L, [PMIN QMIN QMAX]) = [-20 -10 0; -28 0 7];
	t_is( gen, gen1, 8, [t ' - gen'] );
	
	gencost1 = gencost0(:, 1:10);
	gencost1(G, NCOST:(NCOST+6)) = [2 0 0 25 250 0 0; 3 -3 -150 0 0 26 1300; 2 0 0 27 2700 0 0];
	gencost1(L, NCOST:(NCOST+6)) = [[2 -20 -2000 0 0; 2 -28 -280 0 0] zeros(2,2)];
	t_is( gencost, gencost1, 8, [t ' - gencost'] );
	
	t = 'P offers & P bids, lim.P.max_offer/min_bid';
	bids.P.qty = [20; 28];
	bids.P.prc = [100; 10];
	lim.P.min_bid = 50;
	[gen, gencost] = off2case(gen0, gencost0, offers, bids, lim);
	
	gen1 = gen0;
	gen1(G(1:2), PMAX) = offers.P.qty(1:2, :);
	gen1([G(3); L(2)], GEN_STATUS) = 0;
	gen1(L(1), [PMIN QMIN QMAX]) = [-20 -10 0];
	t_is( gen, gen1, 8, [t ' - gen'] );
	
	gencost1 = gencost0;
	gencost1(G(1:2), NCOST:(NCOST+8)) = [[2 0 0 25 250; 2 0 0 26 1300] zeros(2,4)];
	gencost1(L(1), NCOST:(NCOST+8)) = [2 -20 -2000 0 0 0 0 0 0];
	t_is( gencost, gencost1, 8, [t ' - gencost'] );
	
	t = 'P offers & P bids, lim.P.max_offer/min_bid, multi-block';
	offers.P.qty = [10 40; 20 30; 25 25];
	offers.P.prc = [10 100; 25 65; 50 90];
	bids.P.qty = [20 10; 12 18];
	bids.P.prc = [100 60; 70 10];
	[gen, gencost] = off2case(gen0, gencost0, offers, bids, lim);
	
	gen1 = gen0;
	gen1(G, PMAX) = [10; 50; 25];
	gen1(L, [PMIN QMIN QMAX]) = [-30 -15 0; -12 0 3];
	t_is( gen, gen1, 8, [t ' - gen'] );
	
	gencost1 = gencost0(:, 1:10);
	gencost1(G, NCOST:(NCOST+6)) = [2 0 0 10 100 0 0; 3 0 0 20 500 50 2450; 2 0 0 25 1250 0 0];
	gencost1(L, NCOST:(NCOST+6)) = [3 -30 -2600 -20 -2000 0 0; 2 -12 -840 0 0 0 0];
	t_is( gencost, gencost1, 8, [t ' - gencost'] );
	
	%%-----  reactive  -----
	%% generator cost data
	%	1	startup	shutdown	n	x1	y1	...	xn	yn
	%	2	startup	shutdown	n	c(n-1)	...	c0
	gencost0 = [
		1	0	0	4	0	0	12	240	36	1200	60	2400;
		1	100	0	4	0	0	12	240	36	1200	60	2400;
		1	0	0	4	-30	0	-20	1000	-10	2000	0	3000;
		1	0	0	4	0	0	12	240	36	1200	60	2400;
		1	0	50	4	-30	0	-20	1000	-10	2000	0	3000;
		1	0	0	4	-15	-150	0	0	30	150	60	450;
		1	100	0	2	0	0	0	0	0	0	0	0;
		1	0	0	3	-20	-15	-10	-10	0	0	0	0;
		1	0	0	3	0	0	40	80	60	180	0	0;
		1	0	50	2	0	0	0	0	0	0	0	0;
	];
	
	t = 'PQ offers only';
	offers.P.qty = [25; 26; 27];
	offers.P.prc = [10; 50; 100];
	offers.Q.qty = [10; 20; 30];
	offers.Q.prc = [10; 5; 1];
	[gen, gencost] = off2case(gen0, gencost0, offers);
	
	gen1 = gen0;
	gen1(G, PMAX) = offers.P.qty;
	gen1(G, QMAX) = offers.Q.qty;
	gen1(G, QMIN) = 0;
	gen1(L, GEN_STATUS) = 0;
	t_is( gen, gen1, 8, [t ' - gen'] );
	
	gencost1 = gencost0;
	gencost1(G, NCOST:(NCOST+8))     = [[2 0 0 25 250; 2 0 0 26 1300; 2 0 0 27 2700] zeros(3,4)];
	gencost1(G+nGL, NCOST:(NCOST+8)) = [[2 0 0 10 100; 2 0 0 20 100; 2 0 0 30 30] zeros(3,4)];
	
	t_is( gencost, gencost1, 8, [t ' - gencost'] );
	
	t = 'PQ offers & PQ bids, lim.P/Q.max_offer/min_bid, multi-block';
	offers.P.qty = [10 40; 20 30; 25 25];
	offers.P.prc = [10 100; 25 65; 50 90];
	bids.P.qty = [20 10; 12 18];
	bids.P.prc = [100 60; 70 10];
	offers.Q.qty = [5 5; 10 10; 15 15];
	offers.Q.prc = [10 20; 5 60; 1 10];
	bids.Q.qty = [15; 10; 15; 15; 0];
	bids.Q.prc = [-10; 0; 5; -20; 10];
	lim.Q.max_offer = 50;
	lim.Q.min_bid = -15;
	[gen, gencost] = off2case(gen0, gencost0, offers, bids, lim);
	
	gen1 = gen0;
	gen1(:, [GEN_STATUS PMIN PMAX QMIN QMAX]) = [ ...
		1	10	10	-15	10;
		1	12	50	-10	10;
		1	-30	0	-15	0;
		1	12	25	0	30;
		0	-30	0	0	7.5	];
	t_is( gen, gen1, 8, [t ' - gen'] );
	
	gencost1 = gencost0(:, 1:12);
	gencost1(:, NCOST:(NCOST+8)) = [ ...
		2	0	0	10	100	0	0	0	0;
		3	0	0	20	500	50	2450	0	0;
		3	-30	-2600	-20	-2000	0	0	0	0;
		2	0	0	25	1250	0	0	0	0;
		4	-30	0	-20	1000	-10	2000	0	3000;
		4	-15	150	0	0	5	50	10	150;
		3	-10	0	0	0	10	50	0	0;
		2	-15	-75	0	0	0	0	0	0;
		3	0	0	15	15	30	165	0	0;
		2	0	0	0	0	0	0	0	0	];
	t_is( gencost, gencost1, 8, [t ' - gencost'] );
	
	t = 'PQ offers & PQ bids, for gen, no P, no shutdown';
	gen2 = gen0;
	gen2(1, PMIN) = 0;
	offers.P.qty = [0 40; 20 30; 25 25];
	[gen, gencost] = off2case(gen2, gencost0, offers, bids, lim);
	
	gen1(1, [PMIN PMAX QMIN QMAX]) = [ 0 0 -15 10 ];
	t_is( gen, gen1, 8, [t ' - gen'] );
	
	gencost1(1, NCOST:(NCOST+8)) = gencost0(1, NCOST:(NCOST+8));
	t_is( gencost, gencost1, 8, [t ' - gencost'] );
	
	t = 'PQ offers & PQ bids, for gen, no Q, no shutdown';
	offers.P.qty = [10 40; 20 30; 25 25];
	offers.Q.qty = [5 5; 0 10; 15 15];
	bids.Q.qty = [15; 0; 15; 15; 0];
	[gen, gencost] = off2case(gen0, gencost0, offers, bids, lim);
	
	gen1(1, [PMIN PMAX QMIN QMAX]) = [ 10 10 -15 10 ];	%% restore original
	gen1(2, [PMIN PMAX QMIN QMAX]) = [ 12 50 0 0 ];
	t_is( gen, gen1, 8, [t ' - gen'] );
	
	gencost1([1,2,7], NCOST:(NCOST+8)) = [ ...
		2	0	0	10	100	0	0	0	0;
		3	0	0	20	500	50	2450	0	0;
		2	0	0	0	0	0	0	0	0	];
	t_is( gencost, gencost1, 8, [t ' - gencost'] );
	
	t = 'PQ offers & PQ bids, lim.P/Q.max_offer/min_bid, multi-block';
	offers.P.qty = [10 40; 20 30; 25 25];
	offers.P.prc = [10 100; 25 65; 50 90];
	bids.P.qty = [10 0; 12 18];
	bids.P.prc = [100 60; 70 10];
	offers.Q.qty = [5 5; 10 10; 15 15];
	offers.Q.prc = [10 20; 5 60; 1 10];
	bids.Q.qty = [15; 10; 10; 15; 0];
	bids.Q.prc = [-10; 0; 5; -20; 10];
	lim.Q.max_offer = 50;
	lim.Q.min_bid = -15;
	[gen, gencost] = off2case(gen0, gencost0, offers, bids, lim);
	
	gen1 = gen0;
	gen1(:, [GEN_STATUS PMIN PMAX QMIN QMAX]) = [ ...
		1	10	10	-15	10;
		1	12	50	-10	10;
		1	-10	0	-5	0;
		1	12	25	0	30;
		0	-30	0	0	7.5	];
	t_is( gen, gen1, 8, [t ' - gen'] );
	
	gencost1 = gencost0(:, 1:12);
	gencost1(:, NCOST:(NCOST+8)) = [ ...
		2	0	0	10	100	0	0	0	0;
		3	0	0	20	500	50	2450	0	0;
		2	-10	-1000	0	0	0	0	0	0;
		2	0	0	25	1250	0	0	0	0;
		4	-30	0	-20	1000	-10	2000	0	3000;
		4	-15	150	0	0	5	50	10	150;
		3	-10	0	0	0	10	50	0	0;
		2	-10	-50	0	0	0	0	0	0;
		3	0	0	15	15	30	165	0	0;
		2	0	0	0	0	0	0	0	0	];
	t_is( gencost, gencost1, 8, [t ' - gencost'] );

	t = 'PQ offers & PQ bids, zero Q load w/P bid, shutdown bugfix';
	gen1 = gen0;
	gen1(5, [QG, QMIN, QMAX]) = 0;
	[gen, gencost] = off2case(gen1, gencost0, offers, bids, lim);
	
	gen1(:, [PMIN PMAX QMIN QMAX]) = [ ...
		10	10	-15	10;
		12	50	-10	10;
		-10	0	-5	0;
		12	25	0	30;
		-12	0	0	0	];
	t_is( gen, gen1, 8, [t ' - gen'] );
	
	gencost1 = gencost0(:, 1:12);
	gencost1(:, NCOST:(NCOST+8)) = [ ...
		2	0	0	10	100	0	0	0	0;
		3	0	0	20	500	50	2450	0	0;
		2	-10	-1000	0	0	0	0	0	0;
		2	0	0	25	1250	0	0	0	0;
		2	-12	-840	0	0	0	0	0	0;
		4	-15	150	0	0	5	50	10	150;
		3	-10	0	0	0	10	50	0	0;
		2	-10	-50	0	0	0	0	0	0;
		3	0	0	15	15	30	165	0	0;
		2	0	0	0	0	0	0	0	0	];
	t_is( gencost, gencost1, 8, [t ' - gencost'] );

	t = 'PQ offers & PQ bids, non-zero Q load w/no P bid, shutdown bugfix';
	offers.P.qty = [10 40; 20 30; 25 25];
	offers.P.prc = [10 100; 25 65; 50 90];
	bids.P.qty = [0 10; 12 18];
	bids.P.prc = [100 40; 70 10];
	offers.Q.qty = [5 5; 10 10; 15 15];
	offers.Q.prc = [10 20; 5 60; 1 10];
	bids.Q.qty = [15; 10; 15; 15; 0];
	bids.Q.prc = [-10; 0; 5; -20; 10];
	lim.Q.max_offer = 50;
	lim.Q.min_bid = -15;
	[gen, gencost] = off2case(gen0, gencost0, offers, bids, lim);
	
	gen1 = gen0;
	gen1(:, [GEN_STATUS PMIN PMAX QMIN QMAX]) = [ ...
		1	10	10	-15	10;
		1	12	50	-10	10;
		0	-30	0	-15	0;
		1	12	25	0	30;
		0	-30	0	0	7.5	];
	t_is( gen, gen1, 8, [t ' - gen'] );
	
	gencost1 = gencost0(:, 1:12);
	gencost1(:, NCOST:(NCOST+8)) = [ ...
		2	0	0	10	100	0	0	0	0;
		3	0	0	20	500	50	2450	0	0;
		4	-30	0	-20	1000	-10	2000	0	3000;
		2	0	0	25	1250	0	0	0	0;
		4	-30	0	-20	1000	-10	2000	0	3000;
		4	-15	150	0	0	5	50	10	150;
		3	-10	0	0	0	10	50	0	0;
		3	-20	-15	-10	-10	0	0	0	0;
		3	0	0	15	15	30	165	0	0;
		2	0	0	0	0	0	0	0	0	];
	t_is( gencost, gencost1, 8, [t ' - gencost'] );
end

t_end;
