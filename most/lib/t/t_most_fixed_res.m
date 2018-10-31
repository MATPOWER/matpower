function t_most_fixed_res(quiet)
%T_MOST_FIXED_RES  Tests MOST with fixed reserve requirements.

%   MOST
%   Copyright (c) 2012-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

if nargin < 1
    quiet = 0;
end

t_begin(44, quiet);

if quiet
    verbose = 0;
else
    verbose = 0;
end

casefile = 't_case30_userfcns';
mpopt = mpoption('opf.violation', 1e-6, 'mips.gradtol', 1e-8, ...
        'mips.comptol', 1e-8, 'mips.costtol', 1e-9);
mpopt = mpoption(mpopt, 'out.all', 0, 'verbose', verbose, 'opf.ac.solver', 'MIPS');
mpopt = mpoption(mpopt, 'model', 'DC');
mpopt = mpoption(mpopt, 'most.solver', 'DEFAULT');
if have_fcn('gurobi')
    mpopt = mpoption(mpopt, 'gurobi.method', 1);    %% dual-simplex
end
if have_fcn('cplex')
    mpopt = mpoption(mpopt, 'cplex.qpmethod', 2);   %% dual-simplex
end
% if have_fcn('mosek')
%     sc = mosek_symbcon;
%     mpopt = mpoption(mpopt, 'mosek.lp_alg', sc.MSK_OPTIMIZER_DUAL_SIMPLEX);     %% dual simplex
% end
if have_fcn('linprog')
    if have_fcn('linprog_ds')
        mpopt = mpoption(mpopt, 'linprog.Algorithm', 'dual-simplex');
    else
        mpopt = mpoption(mpopt, 'linprog.Algorithm', 'simplex');
    end
end
% mpopt = mpoption(mpopt, 'verbose', 2);

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[CT_LABEL, CT_PROB, CT_TABLE, CT_TBUS, CT_TGEN, CT_TBRCH, CT_TAREABUS, ...
    CT_TAREAGEN, CT_TAREABRCH, CT_ROW, CT_COL, CT_CHGTYPE, CT_REP, ...
    CT_REL, CT_ADD, CT_NEWVAL, CT_TLOAD, CT_TAREALOAD, CT_LOAD_ALL_PQ, ...
    CT_LOAD_FIX_PQ, CT_LOAD_DIS_PQ, CT_LOAD_ALL_P, CT_LOAD_FIX_P, ...
    CT_LOAD_DIS_P, CT_TGENCOST, CT_TAREAGENCOST, CT_MODCOST_F, ...
    CT_MODCOST_X] = idx_ct;

t = 'runopf_w_res(''t_case30_userfcns'') : ';
r1 = runopf_w_res(casefile, mpopt);
t_is(r1.reserves.R, [25; 15; 0; 0; 20; 0], 6, [t 'R']);
t_is(r1.reserves.prc, [2; 2; 2; 2; 5.33007; 5.33007], 5, [t 'prc']);
t_is(r1.reserves.mu.l, [0; 0; 1; 2; 0; 0.169935], 6, [t 'mu.l']);
t_is(r1.reserves.mu.u, [0.1; 0; 0; 0; 0; 0], 7, [t 'mu.u']);
t_is(r1.reserves.mu.Pmax, [0; 0; 0; 0; 0.330065; 0], 6, [t 'mu.Pmax']);
mpc = loadcase(casefile);
[i2e, mpc.bus, mpc.gen, mpc.branch] = ext2int(mpc.bus, mpc.gen, mpc.branch);

t_is(r1.reserves.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(r1.reserves.qty, mpc.reserves.qty, 12, [t 'qty']);
t_is(r1.reserves.totalcost, 177.5, 4, [t 'totalcost']);

%%-----  set up data for DC run (most)  -----
mpc.gen(:, RAMP_10) = Inf;
mpc.gen(:, RAMP_30) = Inf;

%% reserve and delta offers
xgd_table.colnames = {
	'PositiveActiveReservePrice', ...
		'PositiveActiveReserveQuantity', ...
			'NegativeActiveReservePrice', ...
				'NegativeActiveReserveQuantity', ...
					'PositiveActiveDeltaPrice', ...
						'NegativeActiveDeltaPrice', ...
};
xgd_table.data = [
	1	0	1	0	1	1;
	1	0	1	0	1	1;
	1	0	1	0	1	1;
	1	0	1	0	1	1;
	1	0	1	0	1	1;
	1	0	1	0	1	1;
];

ng = size(mpc.gen, 1);      %% number of gens
xgd = loadxgendata(xgd_table, mpc);
mdi = loadmd(mpc, [], xgd);
mdi.FixedReserves = mpc.reserves;

%%-----  run most_fixed_res  -----
%r1 = rundcopf(mpc);
mdo = most(mdi, mpopt);

%%-----  test it  -----
t = 'success1';
t_ok(r1.success, t);

t = 'success2';
t_ok(mdo.QP.exitflag, t);

t = 'f';
t_is(mdo.results.f, r1.f, 7, t);

t = 'Pg';
t_is(mdo.flow.mpc.gen(:, PG), r1.gen(:, PG), 8, t);

t = 'R';
t_is(mdo.flow.mpc.reserves.R, r1.reserves.R, 8, t);

t = 'prc';
t_is(mdo.flow.mpc.reserves.prc, r1.reserves.prc, 8, t);

t = 'totalcost';
t_is(mdo.flow.mpc.reserves.totalcost, r1.reserves.totalcost, 7, t);

t = 'mu.l';
t_is(mdo.flow.mpc.reserves.mu.l, r1.reserves.mu.l, 8, t);

t = 'mu.u';
t_is(mdo.flow.mpc.reserves.mu.u, r1.reserves.mu.u, 8, t);

t = 'mu.Pmax';
t_is(mdo.flow.mpc.reserves.mu.Pmax, r1.reserves.mu.Pmax, 8, t);


%%-----  try again with 3 periods  -----
nt = 3;
profiles = struct( ...
    'type', 'mpcData', ...
    'table', CT_TLOAD, ...
    'rows', 0, ...
    'col', CT_LOAD_ALL_PQ, ...
    'chgtype', CT_REL, ...
    'values', [1.0; 1.1; 1.2] );
mdi = loadmd(mpc, nt, xgd, [], [], profiles);

for t = 1:nt
    mdi.FixedReserves(t,1,1) = mpc.reserves;
end

%%-----  run most  -----
f = 0;
for tt = 1:nt
    mpc1 = mpc;
    mpc1.bus = scale_load(profiles.values(tt), mpc1.bus);
    r(tt) = runopf_w_res(mpc1, mpopt);
    f = f + r(tt).f;
end

mdo = most(mdi, mpopt);

%%-----  test it  -----
t = 'success2';
t_ok(mdo.QP.exitflag, t);

t = 'f';
t_is(mdo.results.f, f, 4, t);

for tt = 1:nt
    t = 'success1';
    t_ok(r(tt).success, t);

    t = 'Pg';
    t_is(mdo.flow(tt,1,1).mpc.gen(:, PG), r(tt).gen(:, PG), 4, sprintf('(t=%d) : %s', tt, t));
    
    t = 'R';
    t_is(mdo.flow(tt,1,1).mpc.reserves.R, r(tt).reserves.R, 4, sprintf('(t=%d) : %s', tt, t));
    
    t = 'prc';
    t_is(mdo.flow(tt,1,1).mpc.reserves.prc, r(tt).reserves.prc, 5, sprintf('(t=%d) : %s', tt, t));
    
    t = 'totalcost';
    t_is(mdo.flow(tt,1,1).mpc.reserves.totalcost, r(tt).reserves.totalcost, 4, sprintf('(t=%d) : %s', tt, t));
    
    t = 'mu.l';
    t_is(mdo.flow(tt,1,1).mpc.reserves.mu.l, r(tt).reserves.mu.l, 5, sprintf('(t=%d) : %s', tt, t));
    
    t = 'mu.u';
    t_is(mdo.flow(tt,1,1).mpc.reserves.mu.u, r(tt).reserves.mu.u, 6, sprintf('(t=%d) : %s', tt, t));
    
    t = 'mu.Pmax';
    t_is(mdo.flow(tt,1,1).mpc.reserves.mu.Pmax, r(tt).reserves.mu.Pmax, 5, sprintf('(t=%d) : %s', tt, t));
end

t_end;
