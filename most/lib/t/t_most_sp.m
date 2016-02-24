function t_most_sp(quiet)
%T_MOST_SP  Tests of single-period continuous optimizations

%   MOST
%   Copyright (c) 2015-2016 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   $Id$
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

t_begin(156, quiet);

if quiet
    verbose = 0;
else
    verbose = 0;
end
% verbose  = 2;

casefile = 't_case3a_most';
mpopt = mpoption;
mpopt = mpoption(mpopt, 'out.gen', 1);
mpopt = mpoption(mpopt, 'verbose', verbose);
% mpopt = mpoption(mpopt, 'opf.violation', 1e-6, 'mips.gradtol', 1e-8, ...
%         'mips.comptol', 1e-8, 'mips.costtol', 1e-8);
mpopt = mpoption(mpopt, 'model', 'DC');
mpopt = mpoption(mpopt, 'sopf.force_Pc_eq_P0', 0);
mpopt = mpoption(mpopt, 'opf.dc.solver', 'MIPS');
mpopt = mpoption(mpopt, 'most.solver', mpopt.opf.dc.solver);
if ~verbose
    mpopt = mpoption(mpopt, 'out.all', 0);
end
% mpopt = mpoption(mpopt, 'out.all', -1);

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[CT_LABEL, CT_PROB, CT_TABLE, CT_TBUS, CT_TGEN, CT_TBRCH, CT_TAREABUS, ...
    CT_TAREAGEN, CT_TAREABRCH, CT_ROW, CT_COL, CT_CHGTYPE, CT_REP, ...
    CT_REL, CT_ADD, CT_NEWVAL, CT_TLOAD, CT_TAREALOAD, CT_LOAD_ALL_PQ, ...
    CT_LOAD_FIX_PQ, CT_LOAD_DIS_PQ, CT_LOAD_ALL_P, CT_LOAD_FIX_P, ...
    CT_LOAD_DIS_P, CT_TGENCOST, CT_TAREAGENCOST, CT_MODCOST_F, ...
    CT_MODCOST_X] = idx_ct;

%% load base case file
mpc = loadcase(casefile);

%% initial xGenData
xgd_table.colnames = {
    'PositiveActiveReservePrice', ...
            'PositiveActiveReserveQuantity', ...
                    'NegativeActiveReservePrice', ...
                            'NegativeActiveReserveQuantity', ...
                                    'PositiveActiveDeltaPrice', ...
                                            'NegativeActiveDeltaPrice', ...
};
xgd_table.data = [
    1e-8    250     2e-8    250     1e-9    1e-9;
    1e-8    250     2e-8    250     1e-9    1e-9;
    1e-8    600     2e-8    600     1e-9    1e-9;
    1e-8    800     2e-8    800     1e-9    1e-9;
];

%%-----  wind  -----
%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
wind.gen = [
	2	0	0	50	-50	1	100	1	100	0	0	0	0	0	0	0	0	200	200	0	0;
];
%% xGenData
wind.xgd_table.colnames = {
	'InitialPg', ...
		'RampWearCostCoeff', ...
			'PositiveActiveReservePrice', ...
				'PositiveActiveReserveQuantity', ...
					'NegativeActiveReservePrice', ...
						'NegativeActiveReserveQuantity', ...
							'PositiveActiveDeltaPrice', ...
								'NegativeActiveDeltaPrice', ...
									'PositiveLoadFollowReservePrice', ...
										'PositiveLoadFollowReserveQuantity', ...
											'NegativeLoadFollowReservePrice', ...
												'NegativeLoadFollowReserveQuantity', ...
};

wind.xgd_table.data = [
	0	0	1e-8	200	2e-8	200	1e-9	1e-9	1e-6	200	1e-6	200;
];

%%-----  contingencies  -----
%% contingency table
% label probty  type        row column      chgtype newvalue
contab = [
    1   0.20    CT_TGEN     2   GEN_STATUS  CT_REP  0;      %% gen 2 at bus 1
    2   0.05    CT_TBRCH    2   BR_STATUS   CT_REP  0;      %% line 1-3
];
pp = [1-sum(contab(:,2)); contab(1,2); contab(2,2)];

xgd = loadxgendata(xgd_table, mpc);

mpc0 = mpc;
xgd0 = xgd;


%%-----  economic dispatch (no network)  -----
if verbose
    fprintf('\n\n');
    fprintf('--------------------------------------------\n');
    fprintf('-----  economic dispatch (no network)  -----\n');
    fprintf('--------------------------------------------\n');
end
t = 'economic dispatch (no network) : runopf ';
mpc.branch(:, RATE_A) = 0;
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.f, -443250, 5, [t 'f']);
t_is(r.gen(:, PG), [150; 150; 150; -450], 7, [t 'Pg']);
t_is(r.bus(:, LAM_P), [30; 30; 30], 7, [t 'lam P']);

%% most
t = 'economic dispatch (no network) : most   ';
mpc = mpc0;
xgd = xgd0;
mpopt = mpoption(mpopt, 'most.dc_model', 0);
mdi = loadmd(mpc, [], xgd);
mdo = most(mdi, mpopt);
rr = mdo.flow(1,1,1).mpc;
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -443250, 5, [t 'f']);
t_is(rr.gen(:, PG), [150; 150; 150; -450], 7, [t 'Pg']);
t_is(rr.bus(:, LAM_P), [30; 30; 30], 7, [t 'lam P']);

%%-----  DC OPF  -----
if verbose
    fprintf('\n\n');
    fprintf('--------------------------------------------\n');
    fprintf('-----              DC OPF              -----\n');
    fprintf('--------------------------------------------\n');
end
t = 'DC OPF : runopf ';
mpc = mpc0;
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.f, -443115, 5, [t 'f']);
t_is(r.gen(:, PG), [135; 135; 180; -450], 7, [t 'Pg']);
t_is(r.bus(:, LAM_P), [27; 36; 45], 7, [t 'lam P']);
t_is(r.branch(:, MU_SF) + r.branch(:, MU_ST), [0; 27; 0], 7, [t 'mu flow']);

%% most
t = 'DC OPF : most   ';
mpc = mpc0;
mpopt = mpoption(mpopt, 'most.dc_model', 1);
mdi = loadmd(mpc, [], xgd);
mdo = most(mdi, mpopt);
rr = mdo.flow(1,1,1).mpc;
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -443115, 5, [t 'f']);
t_is(rr.gen(:, PG), [135; 135; 180; -450], 6, [t 'Pg']);
t_is(rr.bus(:, LAM_P), [27; 36; 45], 7, [t 'lam P']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 27; 0], 6, [t 'mu flow']);


%%-----  economic dispatch (w/reserves)  -----
if verbose
    fprintf('\n\n');
    fprintf('--------------------------------------------\n');
    fprintf('-----  economic dispatch (w/reserves)  -----\n');
    fprintf('--------------------------------------------\n');
end
t = 'economic dispatch (w/reserves) : runopf ';
mpc = mpc0;
mpc.branch(:, RATE_A) = 0;
r = runopf_w_res(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.f, -442820, 5, [t 'f']);
t_is(r.gen(:, PG), [140; 150; 160; -450], 7, [t 'Pg']);
t_is(r.bus(:, LAM_P), [32; 32; 32], 7, [t 'lam P']);
t_is(r.reserves.R, [60; 50; 40; 0], 7, [t 'R']);
t_is(r.reserves.prc, [5; 5; 5; 0], 7, [t 'reserve prc']);
t_is(r.reserves.mu.Pmax + r.gen(:, MU_PMAX), [4; 2; 0; 0], 7, [t 'reserve muPmax']);

%% most
t = 'economic dispatch (w/reserves) : most   ';
mpc = mpc0;
mpopt = mpoption(mpopt, 'most.dc_model', 0);
mdi = loadmd(mpc, [], xgd);
mdi.IncludeFixedReserves = 1;
mdi.FixedReserves = mpc.reserves;
mdo = most(mdi, mpopt);
rr = mdo.flow(1,1,1).mpc;
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -442820, 5, [t 'f']);
t_is(rr.gen(:, PG), [140; 150; 160; -450], 7, [t 'Pg']);
t_is(rr.bus(:, LAM_P), [32; 32; 32], 7, [t 'lam P']);
t_is(mdo.FixedReserves.R, [60; 50; 40; 0], 7, [t 'R']);
t_is(mdo.FixedReserves.prc, [5; 5; 5; 0], 7, [t 'reserve prc']);
t_is(mdo.FixedReserves.mu.Pmax + rr.gen(:, MU_PMAX), [4; 2; 0; 0], 7, [t 'reserve muPmax']);


%%-----  DC OPF  -----
if verbose
    fprintf('\n\n');
    fprintf('--------------------------------------------\n');
    fprintf('-----        DC OPF (w/reserves)       -----\n');
    fprintf('--------------------------------------------\n');
end
t = 'DC OPF (w/reserves) : runopf ';
mpc = mpc0;
r = runopf_w_res(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.f, -442760, 5, [t 'f']);
t_is(r.gen(:, PG), [130; 140; 180; -450], 7, [t 'Pg']);
t_is(r.bus(:, LAM_P), [30; 36; 42], 7, [t 'lam P']);
t_is(r.branch(:, MU_SF) + r.branch(:, MU_ST), [0; 18; 0], 7, [t 'mu flow']);
t_is(r.reserves.R, [70; 60; 20; 0], 7, [t 'R']);
t_is(r.reserves.prc, [5; 5; 5; 0], 7, [t 'reserve prc']);
t_is(r.reserves.mu.Pmax + r.gen(:, MU_PMAX), [4; 2; 0; 0], 7, [t 'reserve muPmax']);

%% most
t = 'DC OPF (w/reserves) : most   ';
mpc = mpc0;
mpopt = mpoption(mpopt, 'most.dc_model', 1);
mdi = loadmd(mpc, [], xgd);
mdi.IncludeFixedReserves = 1;
mdi.FixedReserves = mpc.reserves;
mdo = most(mdi, mpopt);
rr = mdo.flow(1,1,1).mpc;
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -442760, 5, [t 'f']);
t_is(rr.gen(:, PG), [130; 140; 180; -450], 6, [t 'Pg']);
t_is(rr.bus(:, LAM_P), [30; 36; 42], 6, [t 'lam P']);
t_is(mdo.FixedReserves.R, [70; 60; 20; 0], 6, [t 'R']);
t_is(mdo.FixedReserves.prc, [5; 5; 5; 0], 7, [t 'reserve prc']);
t_is(mdo.FixedReserves.mu.Pmax + rr.gen(:, MU_PMAX), [4; 2; 0; 0], 7, [t 'reserve muPmax']);

%%-----  DC OPF w/contingencies  -----
if verbose
    fprintf('\n\n');
    fprintf('--------------------------------------------\n');
    fprintf('-----     DC OPF (w/contingencies)     -----\n');
    fprintf('--------------------------------------------\n');
end
t = 'DC OPF (w/contingencies) : runopf ';
ff = [-443115; -439750; -297000];
mpc = mpc0;
r = runopf(mpc, mpopt);
t_is(r.f, ff(1), 5, [t 'f']);
t_is(r.gen(:, PG), [135; 135; 180; -450], 7, [t 'Pg']);
t_is(r.bus(:, LAM_P), [27; 36; 45], 7, [t 'lam P']);
t_is(r.branch(:, MU_SF) + r.branch(:, MU_ST), [0; 27; 0], 7, [t 'mu flow']);

mpc = mpc0;
mpc.gen(2, GEN_STATUS) = 0;
r = runopf(mpc, mpopt);
t_is(r.f, ff(2), 5, [t 'f']);
t_ok(r.success, [t 'success']);
t_is(r.gen(:, PG), [200; 0; 250; -450], 7, [t 'Pg']);
t_is(r.bus(:, LAM_P), [50; 50; 50], 7, [t 'lam P']);
t_is(r.branch(:, MU_SF) + r.branch(:, MU_ST), [0; 0; 0], 7, [t 'mu flow']);

mpc = mpc0;
mpc.branch(2, BR_STATUS) = 0;
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.f, ff(3), 5, [t 'f']);
t_is(r.gen(:, PG), [100; 100; 100; -300], 7, [t 'Pg']);
t_is(r.bus(:, LAM_P), [20; 20; 1000], 7, [t 'lam P']);
t_is(r.branch(:, MU_SF) + r.branch(:, MU_ST), [0; 0; 980], 7, [t 'mu flow']);

% mpopt = mpoption(mpopt, 'mips.comptol', 1e-8);

t = 'DC OPF (w/contingencies) : c3sopf ';
mpc = mpc0;
r = c3sopf(mpc, xgd_table.data, contab, mpopt);
% r
% r.opf_results
% r.opf_results.f
t_ok(r.success, [t 'success']);
t_is(r.opf_results.f, sum(pp.*ff), 5, [t 'f']);
rr = r.base;
t_is(rr.gen(:, PG), [135; 135; 180; -450], 7, [t 'Pg base']);
t_is(rr.bus(:, LAM_P), [27; 36; 45]*pp(1), 7, [t 'lam P base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 27; 0]*pp(1), 7, [t 'mu flow base']);
rr = r.cont(1);
t_is(rr.gen(:, PG), [200; 0; 250; -450], 7, [t 'Pg 1']);
t_is(rr.bus(:, LAM_P), [50; 50; 50]*pp(2), 7, [t 'lam P 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0]*pp(2), 7, [t 'mu flow 1']);
rr = r.cont(2);
t_is(rr.gen(:, PG), [100; 100; 100; -300], 6, [t 'Pg 2']);
t_is(rr.bus(:, LAM_P), [20; 20; 1000]*pp(3), 7, [t 'lam P 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 980]*pp(3), 7, [t 'mu flow 2']);

t = 'DC OPF (w/contingencies) : most ';
mdi = loadmd(mpc, [], xgd, [], contab);
mdo = most(mdi, mpopt);
% mdo
% mdo.QP
% mdo.QP.f
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -435136.25, 5, [t 'f']);
rr = mdo.flow(1,1,1).mpc;
t_is(rr.gen(:, PG), [135; 135; 180; -450], 7, [t 'Pg base']);
t_is(rr.bus(:, LAM_P), [27; 36; 45]*pp(1), 7, [t 'lam P base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 27; 0]*pp(1), 7, [t 'mu flow base']);
rr = mdo.flow(1,1,2).mpc;
t_is(rr.gen(:, PG), [200; 0; 250; -450], 7, [t 'Pg 1']);
t_is(rr.bus(:, LAM_P), [50; 50; 50]*pp(2), 7, [t 'lam P 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0]*pp(2), 7, [t 'mu flow 1']);
rr = mdo.flow(1,1,3).mpc;
t_is(rr.gen(:, PG), [100; 100; 100; -300], 5, [t 'Pg 2']);
t_is(rr.bus(:, LAM_P), [20; 20; 1000]*pp(3), 7, [t 'lam P 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 980]*pp(3), 7, [t 'mu flow 2']);

% mpopt.verbose = 2;
% mpopt.out.all = -1;

t = 'Secure DC OPF (w/cont,res,ramp) : c3sopf ';
xgd_table.colnames = {
    'PositiveActiveReservePrice', ...
            'PositiveActiveReserveQuantity', ...
                    'NegativeActiveReservePrice', ...
                            'NegativeActiveReserveQuantity', ...
                                    'PositiveActiveDeltaPrice', ...
                                            'NegativeActiveDeltaPrice', ...
};
xgd_table.data = [
    5       250     10      0       1e-9    1e-9;
    1e-8    100     2e-8    0       1e-9    1e-9;
    1.5     600     3       0       1e-9    1e-9;
    1e-8    800     2e-8    0       1e-9    1e-9;
];
xgd = loadxgendata(xgd_table, mpc);
% xgd.PositiveActiveReserveQuantity(2) = 100;
% xgd.PositiveActiveReserveQuantity(4) = 500;
% xgd.NegativeActiveReserveQuantity(:) = 0;
% xgd.PositiveActiveReservePrice([1;3]) = [5; 1.5];
% xgd.NegativeActiveReservePrice([1;3]) = [10; 3];
r = c3sopf(mpc, xgd_table.data, contab, mpopt);
% r
% r.opf_results
% r.opf_results.f
t_ok(r.success, [t 'success']);
t_is(r.opf_results.f, -434638.5625, 5, [t 'f']);
rr = r.base;
% rr.gen(:, PG)
% rr.bus(:, LAM_P)
% rr.branch(:, MU_SF) + rr.branch(:, MU_ST)
t_is(rr.gen(:, PG), [147.5; 122.5; 180; -450], 7, [t 'Pg base']);
t_is(rr.bus(:, LAM_P), [18.8; 27; 35.2], 7, [t 'lam P base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 24.6; 0], 7, [t 'mu flow base']);
rr = r.cont(1);
t_is(rr.gen(:, PG), [181.25; 0; 268.75; -450], 5, [t 'Pg 1']);
t_is(rr.bus(:, LAM_P), [12.25; 12.25; 12.25], 7, [t 'lam P 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 7, [t 'mu flow 1']);
rr = r.cont(2);
t_is(rr.gen(:, PG), [147.5; 22.5; 130; -300], 5, [t 'Pg 2']);
t_is(rr.bus(:, LAM_P), [-0.2; -0.2; 50], 7, [t 'lam P 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 50.2], 7, [t 'mu flow 2']);

t = 'Secure DC OPF (w/cont,res,ramp) : most ';
mdi = loadmd(mpc, [], xgd, [], contab);
mdo = most(mdi, mpopt);
% mdo
% mdo.QP
% mdo.QP.f
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -434638.5625, 4, [t 'f']);
rr = mdo.flow(1,1,1).mpc;
t_is(rr.gen(:, PG), [147.5; 122.5; 180; -450], 6, [t 'Pg base']);
t_is(rr.bus(:, LAM_P), [18.8; 27; 35.2], 7, [t 'lam P base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 24.6; 0], 6, [t 'mu flow base']);
rr = mdo.flow(1,1,2).mpc;
t_is(rr.gen(:, PG), [181.25; 0; 268.75; -450], 5, [t 'Pg 1']);
t_is(rr.bus(:, LAM_P), [12.25; 12.25; 12.25], 7, [t 'lam P 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 7, [t 'mu flow 1']);
rr = mdo.flow(1,1,3).mpc;
t_is(rr.gen(:, PG), [147.5; 22.5; 130; -300], 5, [t 'Pg 2']);
t_is(rr.bus(:, LAM_P), [-0.2; -0.2; 50], 7, [t 'lam P 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 50.2], 7, [t 'mu flow 2']);

t = 'Stochastic DC OPF (w/wind,res) : c3sopf ';
[iwind, mpc, xgd] = addwind(wind, mpc, xgd);
nt = 1;
nj = 3;
% transmat = transmatnormalpersistent(nt, nj);
transmat = {[0.158655253931457; 0.682689492137086; 0.158655253931457]};
profiles = getprofiles(uniformwindprofile(nt, nj), iwind);

pp = transmat{1};
%% contingency table
% label probty  type            row             column          chgtype newvalue
contab0 = [
    1   pp(1)   profiles.table  profiles.rows   profiles.col    profiles.chgtype  profiles.values(1)/profiles.values(2);
    3   pp(3)   profiles.table  profiles.rows   profiles.col    profiles.chgtype  profiles.values(3)/profiles.values(2);
];
mpc0 = mpc;
mpc0.gen(iwind, PMAX) = mpc0.gen(iwind, PMAX) * profiles.values(2);
offers = [xgd_table.data; wind.xgd_table.data(:, [3:8])];
% mpopt.verbose = 2;
% mpopt.out.all = -1;
r = c3sopf(mpc0, offers, contab0, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.opf_results.f, -444525.605052, 5, [t 'f']);
rr = r.cont(1);
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [128.7823267; 141.2176733; 180; -450; 0], 7, [t 'Pg 1']);
t_is(rr.bus(:, LAM_P), [4.4809852; 7.2115891; 9.9421931], 7, [t 'lam P 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 8.1918119; 0], 7, [t 'mu flow 1']);
rr = r.base;
t_is(rr.gen(:, PG), [128.7823267; 135.6088362; 135.6088370; -450; 50], 5, [t 'Pg 2']);
t_is(rr.bus(:, LAM_P), [18.5157455; 18.5157455; 18.5157455], 6, [t 'lam P 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 2']);
rr = r.cont(2);
t_is(rr.gen(:, PG), [128.7823267; 86.9726845; 134.2449888; -450; 100], 5, [t 'Pg 3']);
t_is(rr.bus(:, LAM_P), [2.7597346; 2.7597346; 2.7597346], 6, [t 'lam P 3']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 7, [t 'mu flow 3']);

t = 'Stochastic DC OPF (w/wind,res) : most ';
mdi = loadmd(mpc, transmat, xgd, [], [], profiles);
mdo = most(mdi, mpopt);
% mdo
% mdo.QP
% mdo.QP.f
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -444525.605052, 5, [t 'f']);
rr = mdo.flow(1,1,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [128.7823267; 141.2176733; 180; -450; 0], 7, [t 'Pg 1']);
t_is(rr.bus(:, LAM_P), [4.4809852; 7.2115891; 9.9421931], 7, [t 'lam P 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 8.1918119; 0], 7, [t 'mu flow 1']);
rr = mdo.flow(1,2,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [128.7823267; 135.6088362; 135.6088370; -450; 50], 5, [t 'Pg 2']);
t_is(rr.bus(:, LAM_P), [18.5157455; 18.5157455; 18.5157455], 6, [t 'lam P 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 2']);
rr = mdo.flow(1,3,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [128.7823267; 86.9726845; 134.2449888; -450; 100], 5, [t 'Pg 3']);
t_is(rr.bus(:, LAM_P), [2.7597346; 2.7597346; 2.7597346], 6, [t 'lam P 3']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 7, [t 'mu flow 3']);
% keyboard;

t = 'Secure Stochastic DC OPF (w/wind,cont,res,ramp) : most ';
mdi = loadmd(mpc, transmat, xgd, [], contab, profiles);
mdo = most(mdi, mpopt);
% mdo
% mdo.QP
% mdo.QP.f
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -436240.91814, 3, [t 'f']);

rr = mdo.flow(1,1,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [145.2222662; 124.7777338; 180; -450; 0], 4, [t 'Pg 1 base']);
t_is(rr.bus(:, LAM_P), [3.0669664; 4.2836919; 5.5004173], 3, [t 'lam P 1 base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 3.6501764; 0], 2, [t 'mu flow 1 base']);

rr = mdo.flow(1,1,2).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [156.25; 0; 293.75; -450; 0], 4, [t 'Pg 1 1']);
t_is(rr.bus(:, LAM_P), [3.3641992; 3.3641992; 3.3641992], 7, [t 'lam P 1 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 1 1']);

rr = mdo.flow(1,1,3).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [145.2222662; 42.9397379; 111.8379960; -300; 0], 4, [t 'Pg 1 2']);
t_is(rr.bus(:, LAM_P), [0.0681262; 0.0681262; 7.9327627], 6, [t 'lam P 1 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 7.8646365], 6, [t 'mu flow 1 2']);

rr = mdo.flow(1,2,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [145.2222662; 124.7777338; 130; -450; 50], 3, [t 'Pg 2 base']);
t_is(rr.bus(:, LAM_P), [12.8886634; 13.3124451; 13.7362268], 3, [t 'lam P 2 base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 1.2713451; 0], 2, [t 'mu flow 2 base']);

rr = mdo.flow(1,2,2).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [156.25; 0; 243.75; -450; 50], 3, [t 'Pg 2 1']);
t_is(rr.bus(:, LAM_P), [6.6562225; 6.6562225; 6.6562225], 6, [t 'lam P 2 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 2 1']);

rr = mdo.flow(1,2,3).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [145.2222662; 24.7777338; 111.8379960; -300; 18.1620040], 4, [t 'Pg 2 2']);
t_is(rr.bus(:, LAM_P), [0; 0; 34.1344746], 5, [t 'lam P 2 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 34.1344746], 6, [t 'mu flow 2 2']);

rr = mdo.flow(1,3,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [145.2222662; 92.9397379; 111.8379960; -450; 100], 4, [t 'Pg 3 base']);
t_is(rr.bus(:, LAM_P), [2.2118067; 2.2118067; 2.2118067], 6, [t 'lam P 3 base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 5, [t 'mu flow 3 base']);

rr = mdo.flow(1,3,2).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [156.25; 0; 193.75; -450; 100], 4, [t 'Pg 3 1']);
t_is(rr.bus(:, LAM_P), [1.2295782; 1.2295782; 1.2295782], 6, [t 'lam P 3 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 3 1']);

rr = mdo.flow(1,3,3).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [145.2222662; 24.7777338; 111.8379960; -300; 18.1620040], 4, [t 'Pg 3 2']);
t_is(rr.bus(:, LAM_P), [0; 0; 7.9327628], 5, [t 'lam P 3 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 7.9327627], 6, [t 'mu flow 3 2']);
% keyboard;

t_end;
