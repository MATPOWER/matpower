function t_opf_softlims(quiet)
%T_OPF_SOFTLIMS  Tests for userfcn callbacks (softlims) w/OPF.
%   Includes high-level tests of soft limits implementations.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

casefile = 'case9';
if quiet
    verbose = 0;
else
    verbose = 0;
end

t_begin(59+37*4, quiet);

%% define constants
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% set options
mpopt = mpoption('opf.dc.solver', 'MIPS', 'opf.ac.solver', 'MIPS');
%mpopt = mpoption('opf.dc.solver', 'GLPK', 'opf.ac.solver', 'FMINCON');
%mpopt = mpoption('opf.dc.solver', 'GUROBI', 'opf.ac.solver', 'IPOPT');
%mpopt = mpoption('opf.dc.solver', 'CPLEX', 'opf.ac.solver', 'KNITRO');
%mpopt = mpoption(mpopt, 'fmincon.tol_x', 1e-9, 'fmincon.tol_f', 1e-9);
%mpopt = mpoption(mpopt, 'ipopt.opts.tol', 1e-10, 'ipopt.opts.acceptable_tol', 1e-10);
%mpopt = mpoption(mpopt, 'knitro.tol_x', 1e-9, 'knitro.tol_f', 1e-9);
mpopt = mpoption(mpopt, 'opf.violation', 1e-6, 'mips.gradtol', 1e-8, ...
        'mips.comptol', 1e-8, 'mips.costtol', 1e-9);
if verbose <= 1
    mpopt = mpoption(mpopt, 'out.all', 0);
end
mpopt = mpoption(mpopt, 'verbose', verbose);

%% load and modify case file
mpc = loadcase(casefile);
nl = size(mpc.branch, 1);
mpc.gencost(:, NCOST) = 2;
mpc.gencost(:, COST) = [];
mpc.gencost(:, COST) = [50; 40; 25];
mpc.gencost(:, COST+1) = 0;
mpc.branch(2:6, RATE_A) = [0; 120; 0; 120; 0];
mpc.branch(3, RATE_A) = 120;
mpc.branch(5, RATE_A) = 120;
mpc.branch = [mpc.branch(1:4, :); mpc.branch(4, :);  mpc.branch(5:end, :)];
mpc.branch(5, BR_STATUS) = 0;

nl = size(mpc.branch, 1);   %% number of branches

%% create soft limit inputs
mpc.softlims.idx = (2:nl)';
mpc.softlims.cost = 100 * ones(nl-1, 1);
mpc0 = mpc;     %% save to initialize later runs

t = 'DC - hard limits : ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
r = rundcopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(~isfield(r.softlims, 'overload'), [t 'no softlims.overload']);
t_ok(~isfield(r.softlims, 'ovl_cost'), [t 'no softlims.ovl_cost']);
t_is(r.f, 9126.87, 4, [t 'f']);
t_is(r.gen(:, PG), [12.687; 62.3130; 240], 4, [t 'Pg']);
t_is(r.branch(:, PF), [12.687; -30; -120; 240; 0; 120; 20; -62.3130; 82.3130; -42.687], 4, [t 'Pf']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 35.6504; 0; 0; 7.9756; 0; 0; 0; 0], 4, [t 'mu Pf']);

t = 'DC - soft limits (satisfied) : ';
mpc = toggle_softlims(mpc, 'on');
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = rundcopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims, 'overload'), [t 'softlims.overload exists']);
t_ok(isfield(r.softlims, 'ovl_cost'), [t 'softlims.ovl_cost exists']);
t_is(r.f, 9126.87, 4, [t 'f']);
t_is(r.gen(:, PG), [12.687; 62.3130; 240], 4, [t 'Pg']);
t_is(r.branch(:, PF), [12.687; -30; -120; 240; 0; 120; 20; -62.3130; 82.3130; -42.687], 4, [t 'Pf']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 35.6504; 0; 0; 7.9756; 0; 0; 0; 0], 4, [t 'mu Pf']);
t_is(r.softlims.overload, zeros(10, 1), 12, [t 'softlims.overload']);
t_is(r.softlims.ovl_cost, zeros(10, 1), 12, [t 'softlims.ovl_cost']);
t_is(r.order.branch.status.on(r.order.int.softlims.idx), [3; 6; 8; 9; 10], 12, [t 'mu Pf']);

t = 'savecase(fname, mpc) : ';
fn = sprintf('softlims_savecase_test_%d', fix(1e8*rand));
savecase(fn, mpc);
mpc1 = loadcase(fn);
delete([fn '.m']);
t_ok(isfield(mpc1, 'softlims'), [t 'mpc.softlims']);
t_ok(isfield(mpc1.softlims, 'idx'), [t 'mpc.softlims.idx']);
t_is(mpc1.softlims.idx, mpc.softlims.idx, 5, [t 'mpc.softlims.idx']);
t_ok(isfield(mpc1.softlims, 'cost'), [t 'mpc.softlims.cost']);
t_is(mpc1.softlims.cost, mpc.softlims.cost, 5, [t 'mpc.softlims.cost']);

t = 'savecase(fname, results) + gentype/genfuel: ';
r.genfuel = {'ng'; 'coal'; 'hydro'};
r.gentype = {'GT'; 'ST'; 'HY'};
fn = sprintf('softlims_savecase_test_%d', fix(1e8*rand));
savecase(fn, r);
mpc1 = loadcase(fn);
delete([fn '.m']);
t_ok(isfield(mpc1, 'softlims'), [t 'results.softlims']);
t_ok(isfield(mpc1.softlims, 'idx'), [t 'results.softlims.idx']);
t_is(mpc1.softlims.idx, r.softlims.idx, 5, [t 'results.softlims.idx']);
t_ok(isfield(mpc1.softlims, 'cost'), [t 'results.softlims.cost']);
t_is(mpc1.softlims.cost, r.softlims.cost, 5, [t 'results.softlims.cost']);
t_ok(isfield(mpc1.softlims, 'overload'), [t 'results.softlims.overload']);
t_is(mpc1.softlims.overload, r.softlims.overload, 5, [t 'results.softlims.overload']);
t_ok(isfield(mpc1.softlims, 'ovl_cost'), [t 'results.softlims.ovl_cost']);
t_is(mpc1.softlims.ovl_cost, r.softlims.ovl_cost, 4, [t 'results.softlims.ovl_cost']);
t_ok(isfield(mpc1, 'gentype'), [t 'results.gentype']);
t_ok(isfield(mpc1, 'genfuel'), [t 'results.genfuel']);
if isfield(mpc1, 'gentype')
    t_ok(isequal(mpc1.gentype, r.gentype), [t 'results.gentype']);
else
    t_ok(0, [t 'results.gentype']);
end
if isfield(mpc1, 'genfuel')
    t_ok(isequal(mpc1.genfuel, r.genfuel), [t 'results.genfuel']);
else
    t_ok(0, [t 'results.genfuel']);
end

t = 'DC - soft limits (violated) : ';
mpc = rmfield(mpc, 'softlims');
mpc.softlims.cost = 20;
r = rundcopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims, 'overload'), [t 'softlims.overload exists']);
t_ok(isfield(r.softlims, 'ovl_cost'), [t 'softlims.ovl_cost exists']);
t_is(r.f, 9106.5059, 4, [t 'f']);
t_is(r.gen(:, PG), [10; 63.6988; 241.3012], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -31.3012; -121.3012; 241.3012; 0; 120; 20; -63.6988; 83.6988; -41.3012], 4, [t 'Pf']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 20; 0; 0; 13.2992; 0; 0; 0; 0], 4, [t 'mu Pf']);
t_is(r.softlims.overload, [0; 0; 1.3011811; 0; 0; 0; 0; 0; 0; 0], 6, [t 'softlims.overload']);
t_is(r.softlims.ovl_cost, [0; 0; 26.023622; 0; 0; 0; 0; 0; 0; 0], 6, [t 'softlims.ovl_cost']);

t = 'savecase(fname, mpc) : ';
fn = sprintf('softlims_savecase_test_%d', fix(1e8*rand));
savecase(fn, mpc);
mpc1 = loadcase(fn);
delete([fn '.m']);
t_ok(isfield(mpc1, 'softlims'), [t 'mpc.softlims']);
t_ok(~isfield(mpc1.softlims, 'idx'), [t 'mpc.softlims.idx']);
t_ok(isfield(mpc1.softlims, 'cost'), [t 'mpc.softlims.cost']);
t_is(mpc1.softlims.cost, mpc.softlims.cost, 5, [t 'mpc.softlims.cost']);

t = 'savecase(fname, results) : ';
fn = sprintf('softlims_savecase_test_%d', fix(1e8*rand));
savecase(fn, r);
mpc1 = loadcase(fn);
delete([fn '.m']);
t_ok(isfield(mpc1, 'softlims'), [t 'results.softlims']);
t_ok(isfield(mpc1.softlims, 'idx'), [t 'results.softlims.idx']);
t_ok(isempty(mpc1.softlims.idx), [t 'results.softlims.idx']);
t_ok(isfield(mpc1.softlims, 'cost'), [t 'results.softlims.cost']);
t_is(mpc1.softlims.cost, r.softlims.cost, 5, [t 'results.softlims.cost']);
t_ok(isfield(mpc1.softlims, 'overload'), [t 'results.softlims.overload']);
t_is(mpc1.softlims.overload, r.softlims.overload, 5, [t 'results.softlims.overload']);
t_ok(isfield(mpc1.softlims, 'ovl_cost'), [t 'results.softlims.ovl_cost']);
t_is(mpc1.softlims.ovl_cost, r.softlims.ovl_cost, 4, [t 'results.softlims.ovl_cost']);

%%-----  AC OPF (opf.flow_lim = 'S' -----
mpc = mpc0;
t = 'AC (flow_lim ''S'') - hard limits : ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(~isfield(r.softlims, 'overload'), [t 'no softlims.overload']);
t_ok(~isfield(r.softlims, 'ovl_cost'), [t 'no softlims.ovl_cost']);
t_is(r.f, 9521.385584, 4, [t 'f']);
t_is(r.gen(:, PG), [17.312141; 66.456235; 239.901164], 4, [t 'Pg']);
t_is(r.branch(:, PF), [17.312141; -25.063413; -115.166831; 239.901164; 0; 119.964612; 18.530059; -66.456235; 84.941080; -42.201990], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-17.312141; 25.166831; 119.936552; -239.901164; 0; -118.530059; -18.484844; 66.456235; -82.798010; 42.375554], 4, [t 'Pt']);
t_is(r.branch(:, QF), [-2.721149; -2.842148; -16.241771; 21.275974; 0; -2.914100; -25.325409; -17.123405; 8.944864; -17.404939], 4, [t 'Qf']);
t_is(r.branch(:, QT), [2.884399; -13.758229; -3.901717; 6.815818; 0; -9.674591; 8.178541; 19.603112; -32.595061; -0.042250], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 29.242772; 0; 0; 8.561015; 0; 0; 0; 0], 4, [t 'mu Sf+St']);

t = 'AC (flow_lim ''S'') - soft limits (satisfied) : ';
mpc = toggle_softlims(mpc, 'on');
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims, 'overload'), [t 'softlims.overload exists']);
t_ok(isfield(r.softlims, 'ovl_cost'), [t 'softlims.ovl_cost exists']);
t_is(r.f, 9521.385584, 4, [t 'f']);
t_is(r.gen(:, PG), [17.312141; 66.456235; 239.901164], 4, [t 'Pg']);
t_is(r.branch(:, PF), [17.312141; -25.063413; -115.166831; 239.901164; 0; 119.964612; 18.530059; -66.456235; 84.941080; -42.201990], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-17.312141; 25.166831; 119.936552; -239.901164; 0; -118.530059; -18.484844; 66.456235; -82.798010; 42.375554], 4, [t 'Pt']);
t_is(r.branch(:, QF), [-2.721149; -2.842148; -16.241771; 21.275974; 0; -2.914100; -25.325409; -17.123405; 8.944864; -17.404939], 4, [t 'Qf']);
t_is(r.branch(:, QT), [2.884399; -13.758229; -3.901717; 6.815818; 0; -9.674591; 8.178541; 19.603112; -32.595061; -0.042250], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 29.242772; 0; 0; 8.561015; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.overload, zeros(10, 1), 12, [t 'softlims.overload']);
t_is(r.softlims.ovl_cost, zeros(10, 1), 12, [t 'softlims.ovl_cost']);
t_is(r.order.branch.status.on(r.order.int.softlims.idx), [3; 6; 8; 9; 10], 12, [t 'softlims.idx']);

t = 'AC (flow_lim ''S'') - soft limits (violated) : ';
mpc = rmfield(mpc, 'softlims');
mpc.softlims.cost = 20;
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims, 'overload'), [t 'softlims.overload exists']);
t_ok(isfield(r.softlims, 'ovl_cost'), [t 'softlims.ovl_cost exists']);
t_is(r.f, 9486.265634, 4, [t 'f']);
t_is(r.gen(:, PG), [10; 70.325987; 243.758931], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -28.621625; -118.761245; 243.758931; 0; 119.958848; 18.527937; -70.325987; 88.808549; -38.472172], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-10; 28.761245; 123.8; -243.758931; 0; -118.527937; -18.482561; 70.325987; -86.527828; 38.621625], 4, [t 'Pt']);
t_is(r.branch(:, QF), [3.364131; 0.552145; -12.847555; 19.471020; 0; -3.142423; -25.464124; -14.475416; 6.202740; -20.630111], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-3.306136; -17.152445; -6.346356; 9.488779; 0; -9.535876; 8.272676; 17.182536; -29.369889; 2.753991], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 20; 0; 0; 11.575125; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.overload, [0; 0; 3.962643; 0; 0; 0; 0; 0; 0; 0], 6, [t 'softlims.overload']);
t_is(r.softlims.ovl_cost, [0; 0; 79.252860; 0; 0; 0; 0; 0; 0; 0], 4, [t 'softlims.ovl_cost']);

%%-----  AC OPF (opf.flow_lim = 'I' -----
mpopt = mpoption(mpopt, 'opf.flow_lim', 'I');
mpc = mpc0;
t = 'AC (flow_lim ''I'') - hard limits : ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(~isfield(r.softlims, 'overload'), [t 'no softlims.overload']);
t_ok(~isfield(r.softlims, 'ovl_cost'), [t 'no softlims.ovl_cost']);
t_is(r.f, 9178.300537, 4, [t 'f']);
t_is(r.gen(:, PG), [10; 54.432302; 260.040338], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -32.676659; -122.910718; 260.040338; 0; 131.831065; 30.115930; -54.432302; 84.451894; -42.475047], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-10; 32.910718; 128.209273; -260.040338; 0; -130.115930; -30.019591; 54.432302; -82.524953; 42.676659], 4, [t 'Pt']);
t_is(r.branch(:, QF), [27.824349; 14.237742; 1.341441; 7.733221; 0; -4.887428; -29.462755; -4.743706; -7.797128; -31.778587], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-27.408203; -31.341441; -20.438656; 25.326084; 0; -5.537245; 12.540834; 6.294583; -18.221413; 13.170461], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 0; 0; 0; 20.114712; 0; 0; 0; 0], 4, [t 'mu Sf+St']);

t = 'AC (flow_lim ''I'') - soft limits (satisfied) : ';
mpc = toggle_softlims(mpc, 'on');
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims, 'overload'), [t 'softlims.overload exists']);
t_ok(isfield(r.softlims, 'ovl_cost'), [t 'softlims.ovl_cost exists']);
t_is(r.f, 9178.300537, 4, [t 'f']);
t_is(r.gen(:, PG), [10; 54.432302; 260.040338], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -32.676659; -122.910718; 260.040338; 0; 131.831065; 30.115930; -54.432302; 84.451894; -42.475047], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-10; 32.910718; 128.209273; -260.040338; 0; -130.115930; -30.019591; 54.432302; -82.524953; 42.676659], 4, [t 'Pt']);
t_is(r.branch(:, QF), [27.824349; 14.237742; 1.341441; 7.733221; 0; -4.887428; -29.462755; -4.743706; -7.797128; -31.778587], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-27.408203; -31.341441; -20.438656; 25.326084; 0; -5.537245; 12.540834; 6.294583; -18.221413; 13.170461], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 0; 0; 0; 20.114712; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.overload, zeros(10, 1), 12, [t 'softlims.overload']);
t_is(r.softlims.ovl_cost, zeros(10, 1), 12, [t 'softlims.ovl_cost']);
t_is(r.order.branch.status.on(r.order.int.softlims.idx), [3; 6; 8; 9; 10], 12, [t 'softlims.idx']);

t = 'AC (flow_lim ''I'') - soft limits (violated) : ';
mpc = rmfield(mpc, 'softlims');
mpc.softlims.cost = 15;
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims, 'overload'), [t 'softlims.overload exists']);
t_ok(isfield(r.softlims, 'ovl_cost'), [t 'softlims.ovl_cost exists']);
t_is(r.f, 9146.413838, 4, [t 'f']);
t_is(r.gen(:, PG), [10; 44.862996; 270], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -34.955827; -125.208395; 270; 0; 139.280487; 37.364683; -44.862996; 82.094229; -44.742338], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-10; 35.208395; 130.719513; -270; 0; -137.364683; -37.231233; 44.862996; -80.257662; 44.955827], 4, [t 'Pt']);
t_is(r.branch(:, QF), [25.262645; 13.437313; 0.303801; 13.552057; 0; -3.645909; -29.947558; -7.234346; -6.145528; -29.826268], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-24.907470; -30.303801; -18.341454; 21.987363; 0; -5.052442; 13.379873; 8.309624; -20.173732; 11.470157], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 6.260861; 0; 0; 15; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.overload, [0; 0; 0; 0; 0; 6.792934; 0; 0; 0; 0], 6, [t 'softlims.overload']);
t_is(r.softlims.ovl_cost, [0; 0; 0; 0; 0; 101.894009; 0; 0; 0; 0], 4, [t 'softlims.ovl_cost']);

%%-----  AC OPF (opf.flow_lim = '2' -----
mpopt = mpoption(mpopt, 'opf.flow_lim', '2');
mpc = mpc0;
t = 'AC (flow_lim ''2'') - hard limits : ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(~isfield(r.softlims, 'overload'), [t 'no softlims.overload']);
t_ok(~isfield(r.softlims, 'ovl_cost'), [t 'no softlims.ovl_cost']);
t_is(r.f, 9513.613051, 4, [t 'f']);
t_is(r.gen(:, PG), [17.759246; 65.641269; 240], 4, [t 'Pg']);
t_is(r.branch(:, PF), [17.759246; -25.216499; -115.346644; 240; 0; 120; 18.577311; -65.641269; 84.170657; -42.784248], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-17.759246; 25.346644; 120; -240; 0; -118.577311; -18.529387; 65.641269; -82.215752; 42.975745], 4, [t 'Pt']);
t_is(r.branch(:, QF), [15.646152; 6.577259; -6.087866; 8.021258; 0; -4.520662; -26.623563; -6.687069; -2.629290; -27.085546], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-15.370220; -23.912134; -15.548682; 20.069344; 0; -8.376437; 9.316359; 8.954090; -22.914454; 8.792961], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 29.248033; 0; 0; 8.417115; 0; 0; 0; 0], 4, [t 'mu Sf+St']);

t = 'AC (flow_lim ''2'') - soft limits (satisfied) : ';
mpc = toggle_softlims(mpc, 'on');
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims, 'overload'), [t 'softlims.overload exists']);
t_ok(isfield(r.softlims, 'ovl_cost'), [t 'softlims.ovl_cost exists']);
t_is(r.f, 9513.613051, 4, [t 'f']);
t_is(r.gen(:, PG), [17.759246; 65.641269; 240], 4, [t 'Pg']);
t_is(r.branch(:, PF), [17.759246; -25.216499; -115.346644; 240; 0; 120; 18.577311; -65.641269; 84.170657; -42.784248], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-17.759246; 25.346644; 120; -240; 0; -118.577311; -18.529387; 65.641269; -82.215752; 42.975745], 4, [t 'Pt']);
t_is(r.branch(:, QF), [15.646152; 6.577259; -6.087866; 8.021258; 0; -4.520662; -26.623563; -6.687069; -2.629290; -27.085546], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-15.370220; -23.912134; -15.548682; 20.069344; 0; -8.376437; 9.316359; 8.954090; -22.914454; 8.792961], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 29.248033; 0; 0; 8.417115; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.overload, zeros(10, 1), 12, [t 'softlims.overload']);
t_is(r.softlims.ovl_cost, zeros(10, 1), 12, [t 'softlims.ovl_cost']);
t_is(r.order.branch.status.on(r.order.int.softlims.idx), [3; 6; 8; 9; 10], 12, [t 'softlims.idx']);

t = 'AC (flow_lim ''2'') - soft limits (violated) : ';
mpc = rmfield(mpc, 'softlims');
mpc.softlims.cost = 20;
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims, 'overload'), [t 'softlims.overload exists']);
t_ok(isfield(r.softlims, 'ovl_cost'), [t 'softlims.ovl_cost exists']);
t_is(r.f, 9476.854715, 4, [t 'f']);
t_is(r.gen(:, PG), [10; 69.833910; 244.077740; ], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -28.933383; -119.111239; 244.077740; 0; 120; 18.578292; -69.833910; 88.362736; -38.762267], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-10; 29.111239; 124.077740; -244.077740; 0; -118.578292; -18.528826; 69.833910; -86.237733; 38.933383], 4, [t 'Pt']);
t_is(r.branch(:, QF), [22.204979; 10.325903; -2.434897; 5.980065; 0; -5.159694; -27.240912; -4.815545; -5.105567; -30.250169], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-21.917920; -27.565103; -17.970864; 23.130558; 0; -7.759088; 9.921112; 7.362539; -19.749831; 11.592017], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 20; 0; 0; 11.526783; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.overload, [0; 0; 4.07774; 0; 0; 0; 0; 0; 0; 0], 6, [t 'softlims.overload']);
t_is(r.softlims.ovl_cost, [0; 0; 81.554809; 0; 0; 0; 0; 0; 0; 0], 4, [t 'softlims.ovl_cost']);

%%-----  AC OPF (opf.flow_lim = 'P' -----
mpopt = mpoption(mpopt, 'opf.flow_lim', 'P');
mpc = mpc0;
t = 'AC (flow_lim ''P'') - hard limits : ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(~isfield(r.softlims, 'overload'), [t 'no softlims.overload']);
t_ok(~isfield(r.softlims, 'ovl_cost'), [t 'no softlims.ovl_cost']);
t_is(r.f, 9513.613051, 4, [t 'f']);
t_is(r.gen(:, PG), [17.759246; 65.641269; 240], 4, [t 'Pg']);
t_is(r.branch(:, PF), [17.759246; -25.216499; -115.346644; 240; 0; 120; 18.577311; -65.641269; 84.170657; -42.784248], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-17.759246; 25.346644; 120; -240; 0; -118.577311; -18.529387; 65.641269; -82.215752; 42.975745], 4, [t 'Pt']);
t_is(r.branch(:, QF), [15.646152; 6.577259; -6.087866; 8.021258; 0; -4.520662; -26.623563; -6.687069; -2.629290; -27.085546], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-15.370220; -23.912134; -15.548682; 20.069344; 0; -8.376437; 9.316359; 8.954090; -22.914454; 8.792961], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 29.248033; 0; 0; 8.417115; 0; 0; 0; 0], 4, [t 'mu Sf+St']);

t = 'AC (flow_lim ''P'') - soft limits (satisfied) : ';
mpc = toggle_softlims(mpc, 'on');
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims, 'overload'), [t 'softlims.overload exists']);
t_ok(isfield(r.softlims, 'ovl_cost'), [t 'softlims.ovl_cost exists']);
t_is(r.f, 9513.613051, 4, [t 'f']);
t_is(r.gen(:, PG), [17.759246; 65.641269; 240], 4, [t 'Pg']);
t_is(r.branch(:, PF), [17.759246; -25.216499; -115.346644; 240; 0; 120; 18.577311; -65.641269; 84.170657; -42.784248], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-17.759246; 25.346644; 120; -240; 0; -118.577311; -18.529387; 65.641269; -82.215752; 42.975745], 4, [t 'Pt']);
t_is(r.branch(:, QF), [15.646152; 6.577259; -6.087866; 8.021258; 0; -4.520662; -26.623563; -6.687069; -2.629290; -27.085546], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-15.370220; -23.912134; -15.548682; 20.069344; 0; -8.376437; 9.316359; 8.954090; -22.914454; 8.792961], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 29.248033; 0; 0; 8.417115; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.overload, zeros(10, 1), 12, [t 'softlims.overload']);
t_is(r.softlims.ovl_cost, zeros(10, 1), 12, [t 'softlims.ovl_cost']);
t_is(r.order.branch.status.on(r.order.int.softlims.idx), [3; 6; 8; 9; 10], 12, [t 'softlims.idx']);

t = 'AC (flow_lim ''P'') - soft limits (violated) : ';
mpc = rmfield(mpc, 'softlims');
mpc.softlims.cost = 20;
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims, 'overload'), [t 'softlims.overload exists']);
t_ok(isfield(r.softlims, 'ovl_cost'), [t 'softlims.ovl_cost exists']);
t_is(r.f, 9476.854715, 4, [t 'f']);
t_is(r.gen(:, PG), [10; 69.833910; 244.077740; ], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -28.933383; -119.111239; 244.077740; 0; 120; 18.578292; -69.833910; 88.362736; -38.762267], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-10; 29.111239; 124.077740; -244.077740; 0; -118.578292; -18.528826; 69.833910; -86.237733; 38.933383], 4, [t 'Pt']);
t_is(r.branch(:, QF), [22.204979; 10.325903; -2.434897; 5.980065; 0; -5.159694; -27.240912; -4.815545; -5.105567; -30.250169], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-21.917920; -27.565103; -17.970864; 23.130558; 0; -7.759088; 9.921112; 7.362539; -19.749831; 11.592017], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 20; 0; 0; 11.526783; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.overload, [0; 0; 4.07774; 0; 0; 0; 0; 0; 0; 0], 6, [t 'softlims.overload']);
t_is(r.softlims.ovl_cost, [0; 0; 81.554809; 0; 0; 0; 0; 0; 0; 0], 4, [t 'softlims.ovl_cost']);

t_end;
