function t_runmarket(quiet)
%T_RUNMARKET  Tests for code in runmkt, smartmkt and auction.m.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2005 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

n_tests = 34;

t_begin(n_tests, quiet);

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

if ~have_fcn('smartmarket')
    t_skip(n_tests, 'smartmarket code not available');
elseif ~have_fcn('minopf')
    t_skip(n_tests, 't_runmarket requires MINOPF');
else
	mpc = loadcase('t_auction_case');
	
	mpopt = mpoption('OPF_ALG', 500, 'OUT_ALL_LIM', 1, 'OUT_BRANCH', 0, 'OUT_SYS_SUM', 0, 'OUT_ALL', 0, 'VERBOSE', 1);
	% mpopt = mpoption('OUT_GEN', 1, 'OUT_BRANCH', 0, 'OUT_SYS_SUM', 0);
	
	offers.P.qty = [
		12 24 24; 
		12 24 24; 
		12 24 24; 
		12 24 24; 
		12 24 24; 
		12 24 24; 
	];
	offers.P.prc = [
		20 50 60;
		20 40 70;
		20 42 80;
		20 44 90;
		20 46 75;
		20 48 60;
	];
	bids.P.qty = [
		10 10 10;
		10 10 10;
		10 10 10;
	];
	bids.P.prc = [
		100 70 60;
	% 	100 64.3 20;
	% 	100 30.64545 0;
		100 50 20;
		100 60 50;
	];
	
	offers.Q.qty = [ 60; 60; 60; 60; 60; 60; 0; 0; 0 ];
	offers.Q.prc = [ 0; 0; 0; 0; 0; 3; 0; 0; 0 ];
	bids.Q.qty = [ 15; 15; 15; 15; 15; 15; 15; 12; 7.5 ];
	% bids.Q.prc = [ 0; 0; 0; 0; 0; 0; 0; 83.9056; 0 ];
	bids.Q.prc = [ 0; 0; 0; 0; 0; 0; 0; 20; 0 ];
	
	t = 'marginal Q offer, marginal PQ bid, auction_type = 5';
	mkt = struct(	'auction_type', 5, ...
					't' , [], ...
					'u0', [], ...
					'lim', []	);
	[c, co, cb, f, dispatch, success, et] = runmarket(mpc, offers, bids, mkt, mpopt);
	co5 = co;
	cb5 = cb;
	
	% [ co.P.qty co.P.prc ]
	% [ cb.P.qty cb.P.prc ]
	% [ co.Q.qty co.Q.prc ]
	% [ cb.Q.qty cb.Q.prc ]
	
	G = find( ~isload(c.gen) );   %% real generators
	L = find(  isload(c.gen) );   %% dispatchable loads
	Gbus = c.gen(G,GEN_BUS);
	Lbus = c.gen(L,GEN_BUS);
	
	t_is( co.P.qty, ones(6, 1) * [12 24 0], 2, [t ' : gen P quantities'] );
	t_is( co.P.prc(1,:), 50.1578, 4, [t ' : gen 1 P prices'] );
	t_is( cb.P.qty, [10 10 10; 10 0.196 0; 10 10 0], 2, [t ' : load P quantities'] );
	t_is( cb.P.prc(2,:), 56.9853, 4, [t ' : load 2 P price'] );
	t_is( co.P.prc(:,1), c.bus(Gbus, LAM_P), 8, [t ' : gen P prices'] );
	t_is( cb.P.prc(:,1), c.bus(Lbus, LAM_P), 8, [t ' : load P prices'] );
	
	lao_gap   = offers.P.prc(1,2) - c.bus(Gbus(1), LAM_P);
	fro_gap   = offers.P.prc(6,3) - c.bus(Gbus(6), LAM_P);
	
	t_is( lao_gap, -0.1578, 4, 'P lao_gap');
	t_is( fro_gap, 3.9802, 4, 'P fro_gap');
	
	t_is( co.Q.qty, [4.2722; 11.3723; 14.1472; 22.8939; 36.7886; 12.3375; 0; 0; 0], 2, [t ' : Q offer quantities'] );
	t_is( co.Q.prc, [0;0;0;0;0;3; 0.4861; 2.5367; 1.3763], 4, [t ' : Q offer prices'] );
	t_is( cb.Q.qty, [0;0;0;0;0;0; 15; 4.0785; 5], 2, [t ' : Q bid quantities'] );
	t_is( cb.Q.prc, [0;0;0;0;0;3; 0.4861; 2.5367; 1.3763], 4, [t ' : Q bid prices'] );
	t_is( co.Q.prc, c.bus([Gbus; Lbus], LAM_Q), 8, [t ' : Q offer prices'] );
	t_is( cb.Q.prc, co.Q.prc, 8, [t ' : Q bid prices'] );
	
	t = 'marginal Q offer, marginal PQ bid, auction_type = 0';
	mkt.auction_type = 0;
	[c, co, cb, f, dispatch, success, et] = runmarket(mpc, offers, bids, mkt, mpopt);
	t_is( co.P.qty, co5.P.qty, 8, [t ' : gen P quantities'] );
	t_is( cb.P.qty, cb5.P.qty, 8, [t ' : load P quantities'] );
	t_is( co.P.prc, offers.P.prc, 8, [t ' : gen P prices'] );
	t_is( cb.P.prc, bids.P.prc, 8, [t ' : load P prices'] );
	
	t_is( co.Q.qty, co5.Q.qty, 8, [t ' : gen Q quantities'] );
	t_is( cb.Q.qty, cb5.Q.qty, 8, [t ' : load Q quantities'] );
	t_is( co.Q.prc, offers.Q.prc, 8, [t ' : gen Q prices'] );
	t_is( cb.Q.prc, bids.Q.prc, 8, [t ' : load Q prices'] );
	
	t = 'marginal Q offer, marginal PQ bid, auction_type = 1';
	mkt.auction_type = 1;
	[c, co, cb, f, dispatch, success, et] = runmarket(mpc, offers, bids, mkt, mpopt);
	
	t_is( co.P.qty, co5.P.qty, 8, [t ' : gen P quantities'] );
	t_is( cb.P.qty, cb5.P.qty, 8, [t ' : load P quantities'] );
	
	t_is( co.P.prc(1,:), 50, 4, [t ' : gen 1 P prices'] );
	t_is( cb.P.prc, cb5.P.prc + lao_gap, 4, [t ' : load 2 P price'] );
	t_is( co.P.prc(:,1), c.bus(Gbus, LAM_P) + lao_gap, 8, [t ' : gen P prices'] );
	t_is( cb.P.prc(:,1), c.bus(Lbus, LAM_P) + lao_gap, 8, [t ' : load P prices'] );
	
	t_is( co.Q.qty, co5.Q.qty, 8, [t ' : gen Q quantities'] );
	t_is( cb.Q.qty, cb5.Q.qty, 8, [t ' : load Q quantities'] );
	
	t_is( co.Q.qty, co5.Q.qty, 4, [t ' : Q offer quantities'] );
	t_is( co.Q.prc, co5.Q.prc, 4, [t ' : Q offer prices'] );
	t_is( cb.Q.qty, cb5.Q.qty, 4, [t ' : Q bid quantities'] );
	t_is( cb.Q.prc, cb5.Q.prc, 4, [t ' : Q bid prices'] );
end

t_end;

return;
