function t_most_spuc(quiet)
%T_MOST_SPUC  Tests of single-period unit commitment optimizations

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

solvers = {'CPLEX', 'GLPK', 'GUROBI', 'MOSEK', 'OT'};
fcn = {'cplex', 'glpk', 'gurobi', 'mosek', 'intlinprog'};
% solvers = {'OT'};
% fcn = {'intlinprog'};
% solvers = {'GUROBI'};
% fcn = {'gurobi'};
% solvers = {'GLPK'};
% fcn = {'glpk'};
ntests = 144;
t_begin(ntests*length(solvers), quiet);

if quiet
    verbose = 0;
else
    verbose = 0;
end
% verbose = 2;

casefile = 'ex_case3b';
%solnfile =  't_most_spuc_soln';
%soln = load(solnfile);
mpopt = mpoption;
mpopt = mpoption(mpopt, 'out.gen', 1);
mpopt = mpoption(mpopt, 'verbose', verbose);
% mpopt = mpoption(mpopt, 'opf.violation', 1e-6, 'mips.gradtol', 1e-8, ...
%         'mips.comptol', 1e-8, 'mips.costtol', 1e-8);
mpopt = mpoption(mpopt, 'model', 'DC');
mpopt = mpoption(mpopt, 'sopf.force_Pc_eq_P0', 0);
mpopt = mpoption(mpopt, 'most.price_stage_warn_tol', 10);

%% solver options
if have_fcn('cplex')
    %mpopt = mpoption(mpopt, 'cplex.lpmethod', 0);       %% automatic
    %mpopt = mpoption(mpopt, 'cplex.lpmethod', 1);       %% primal simplex
    mpopt = mpoption(mpopt, 'cplex.lpmethod', 2);       %% dual simplex
    %mpopt = mpoption(mpopt, 'cplex.lpmethod', 3);       %% network simplex
    %mpopt = mpoption(mpopt, 'cplex.lpmethod', 4);       %% barrier
    mpopt = mpoption(mpopt, 'cplex.opts.mip.tolerances.mipgap', 0);
    mpopt = mpoption(mpopt, 'cplex.opts.mip.tolerances.absmipgap', 0);
    mpopt = mpoption(mpopt, 'cplex.opts.threads', 2);
end
if have_fcn('glpk')
    mpopt = mpoption(mpopt, 'glpk.opts.mipgap', 0);
    mpopt = mpoption(mpopt, 'glpk.opts.tolint', 1e-10);
    mpopt = mpoption(mpopt, 'glpk.opts.tolobj', 1e-10);
end
if have_fcn('gurobi')
    %mpopt = mpoption(mpopt, 'gurobi.method', -1);       %% automatic
    %mpopt = mpoption(mpopt, 'gurobi.method', 0);        %% primal simplex
    mpopt = mpoption(mpopt, 'gurobi.method', 1);        %% dual simplex
    %mpopt = mpoption(mpopt, 'gurobi.method', 2);        %% barrier
    mpopt = mpoption(mpopt, 'gurobi.threads', 2);
    mpopt = mpoption(mpopt, 'gurobi.opts.MIPGap', 0);
    mpopt = mpoption(mpopt, 'gurobi.opts.MIPGapAbs', 0);
end
if have_fcn('mosek')
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', 0);         %% automatic
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', 1);         %% interior point
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', 3);         %% primal simplex
    mpopt = mpoption(mpopt, 'mosek.lp_alg', 4);         %% dual simplex
    sc = mosek_symbcon;
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', 5);         %% primal dual simplex
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', 6);         %% automatic simplex
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', 7);         %% network primal simplex
    %mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_MIO_TOL_X', 0);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_IPAR_MIO_NODE_OPTIMIZER', sc.MSK_OPTIMIZER_DUAL_SIMPLEX);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_IPAR_MIO_ROOT_OPTIMIZER', sc.MSK_OPTIMIZER_DUAL_SIMPLEX);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_MIO_TOL_ABS_RELAX_INT', 1e-9);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_MIO_TOL_REL_RELAX_INT', 0);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_MIO_TOL_REL_GAP', 0);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_MIO_TOL_ABS_GAP', 0);
end
if have_fcn('intlinprog')
    %mpopt = mpoption(mpopt, 'linprog.Algorithm', 'interior-point');
    %mpopt = mpoption(mpopt, 'linprog.Algorithm', 'active-set');
    %mpopt = mpoption(mpopt, 'linprog.Algorithm', 'simplex');
    mpopt = mpoption(mpopt, 'linprog.Algorithm', 'dual-simplex');
    %mpopt = mpoption(mpopt, 'intlinprog.RootLPAlgorithm', 'primal-simplex');
    mpopt = mpoption(mpopt, 'intlinprog.RootLPAlgorithm', 'dual-simplex');
    mpopt = mpoption(mpopt, 'intlinprog.TolCon', 1e-9);
    mpopt = mpoption(mpopt, 'intlinprog.TolGapAbs', 0);
    mpopt = mpoption(mpopt, 'intlinprog.TolGapRel', 0);
    mpopt = mpoption(mpopt, 'intlinprog.TolInteger', 1e-6);
end
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
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
[CT_LABEL, CT_PROB, CT_TABLE, CT_TBUS, CT_TGEN, CT_TBRCH, CT_TAREABUS, ...
    CT_TAREAGEN, CT_TAREABRCH, CT_ROW, CT_COL, CT_CHGTYPE, CT_REP, ...
    CT_REL, CT_ADD, CT_NEWVAL, CT_TLOAD, CT_TAREALOAD, CT_LOAD_ALL_PQ, ...
    CT_LOAD_FIX_PQ, CT_LOAD_DIS_PQ, CT_LOAD_ALL_P, CT_LOAD_FIX_P, ...
    CT_LOAD_DIS_P, CT_TGENCOST, CT_TAREAGENCOST, CT_MODCOST_F, ...
    CT_MODCOST_X] = idx_ct;

%% load base case file
mpc = loadcase(casefile);
mpc.gencost(:, STARTUP) = 0;
mpc.gencost(:, SHUTDOWN) = 0;

%%-----  contingencies  -----
contab = loadgenericdata('ex_contab', 'array');
pp = [1-sum(contab(:,2)); contab(1,2); contab(2,2)];

xgd = loadxgendata('ex_xgd', mpc);
[iwind, mpc, xgd] = addwind('ex_wind', mpc, xgd);
mpc.reserves.zones = [mpc.reserves.zones 0];
mpc.gen(4, PMIN) = -499;
mpc0 = mpc;
xgd0 = xgd;

for s = 1:length(solvers)
    if ~have_fcn(fcn{s})     %% check if we have the solver
        t_skip(ntests, sprintf('%s not installed', solvers{s}));
    else
        mpopt = mpoption(mpopt, 'opf.dc.solver', solvers{s});
        mpopt = mpoption(mpopt, 'most.solver', mpopt.opf.dc.solver);

%%-----  economic dispatch (no network)  -----
if verbose
    fprintf('\n\n');
    fprintf('--------------------------------------------\n');
    fprintf('-----  economic dispatch (no network)  -----\n');
    fprintf('--------------------------------------------\n');
end
t = sprintf('%s : economic dispatch (no network) : runopf ', solvers{s});
mpc = mpc0;
mpc.branch(:, RATE_A) = 0;
r = runuopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.f, -488030, 7, [t 'f']);
t_is(r.gen(:, PG), [200; 199; 0; -499; 100], 7, [t 'Pg']);
t_is(r.gen(:, GEN_STATUS), [1; 1; 0; 1; 1], 7, [t 'u']);
t_is(r.bus(:, LAM_P), [30; 30; 30], 7, [t 'lam P']);

%% most
t = sprintf('%s : economic dispatch (no network) : most   ', solvers{s});
mpc = mpc0;
mpc.gen(1, GEN_STATUS) = 0;
mpc.gen(2, GEN_STATUS) = 0;
% xgd = xg
mpopt = mpoption(mpopt, 'most.dc_model', 0);
mdi = loadmd(mpc, [], xgd);
mdo = most(mdi, mpopt);
rr = mdo.flow(1,1,1).mpc;
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -488030, 7, [t 'f']);
t_is(rr.gen(:, PG), [200; 199; 0; -499; 100], 7, [t 'Pg']);
t_is(rr.gen(:, GEN_STATUS), [1; 1; 0; 1; 1], 7, [t 'u']);
% rr.gen(:, GEN_STATUS)
t_is(rr.bus(:, LAM_P), [30; 30; 30], 7, [t 'lam P']);

%%-----  DC OPF  -----
if verbose
    fprintf('\n\n');
    fprintf('--------------------------------------------\n');
    fprintf('-----              DC OPF              -----\n');
    fprintf('--------------------------------------------\n');
end
t = sprintf('%s : DC OPF : runopf ', solvers{s});
mpc = mpc0;
% mpc.gen(4, PMIN) = -500;
r = runuopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.f, -486040, 7, [t 'f']);
t_is(r.gen(:, PG), [200; 0; 199; -499; 100], 7, [t 'Pg']);
t_is(r.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(r.bus(:, LAM_P), [40; 40; 40], 7, [t 'lam P']);
t_is(r.branch(:, MU_SF) + r.branch(:, MU_ST), [0; 0; 0], 7, [t 'mu flow']);

%% most
t = sprintf('%s : DC OPF : most   ', solvers{s});
mpc = mpc0;
% mpc.gen(4, PMIN) = -500;
mpopt = mpoption(mpopt, 'most.dc_model', 1);
mdi = loadmd(mpc, [], xgd);
mdo = most(mdi, mpopt);
rr = mdo.flow(1,1,1).mpc;
if verbose
    printpf(rr, [], mpopt);
end
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -486040, 7, [t 'f']);
t_is(rr.gen(:, PG), [200; 0; 199; -499; 100], 6, [t 'Pg']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
% rr.gen(:, GEN_STATUS)
t_is(rr.bus(:, LAM_P), [40; 40; 40], 7, [t 'lam P']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow']);

%%-----  economic dispatch (w/reserves)  -----
if verbose
    fprintf('\n\n');
    fprintf('--------------------------------------------\n');
    fprintf('-----  economic dispatch (w/reserves)  -----\n');
    fprintf('--------------------------------------------\n');
end
t = sprintf('%s : economic dispatch (w/reserves) : runopf ', solvers{s});
mpc = mpc0;
mpc.branch(:, RATE_A) = 0;
mpc = toggle_reserves(mpc, 'on');
r = runuopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.f, -486802, 5, [t 'f']);
t_is(r.gen(:, PG), [200; 139; 60; -499; 100], 7, [t 'Pg']);
t_is(r.gen(:, GEN_STATUS), [1; 1; 1; 1; 1], 7, [t 'u']);
t_is(r.bus(:, LAM_P), [32; 32; 32], 7, [t 'lam P']);
t_is(r.reserves.R, [0; 61; 89; 0; 0], 7, [t 'R']);
t_is(r.reserves.prc, [5; 5; 5; 0; 0], 7, [t 'reserve prc']);
t_is(r.reserves.mu.Pmax + r.gen(:, MU_PMAX), [7; 2; 0; 0; 32], 7, [t 'reserve muPmax']);

%% most
t = sprintf('%s : economic dispatch (w/reserves) : most   ', solvers{s});
mpc = mpc0;
% mpc.gen(4, PMIN) = -350.8;
mpopt = mpoption(mpopt, 'most.dc_model', 0);
mdi = loadmd(mpc, [], xgd);
mdi.FixedReserves = mpc.reserves;
mdo = most(mdi, mpopt);
rr = mdo.flow(1,1,1).mpc;
if verbose
    printpf(rr, [], mpopt);
end
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -486802, 5, [t 'f']);
t_is(rr.gen(:, PG), [200; 139; 60; -499; 100], 7, [t 'Pg']);
t_is(rr.gen(:, GEN_STATUS), [1; 1; 1; 1; 1], 7, [t 'u']);
% rr.gen(:, GEN_STATUS)
t_is(rr.bus(:, LAM_P), [32; 32; 32], 7, [t 'lam P']);
t_is(rr.reserves.R, [0; 61; 89; 0; 0], 7, [t 'R']);
t_is(rr.reserves.prc, [5; 5; 5; 0; 0], 7, [t 'reserve prc']);
t_is(rr.reserves.mu.Pmax + rr.gen(:, MU_PMAX), [7; 2; 0; 0; 32], 7, [t 'reserve muPmax']);


%%-----  DC OPF  -----
if verbose
    fprintf('\n\n');
    fprintf('--------------------------------------------\n');
    fprintf('-----        DC OPF (w/reserves)       -----\n');
    fprintf('--------------------------------------------\n');
end
t = sprintf('%s : DC OPF (w/reserves) : runopf ', solvers{s});
mpc = mpc0;
r = runopf_w_res(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.f, -485656, 5, [t 'f']);
t_is(r.gen(:, PG), [156; 65; 178; -499; 100], 7, [t 'Pg']);
t_is(r.gen(:, GEN_STATUS), [1; 1; 1; 1; 1], 7, [t 'u']);
t_is(r.bus(:, LAM_P), [29; 40; 51], 7, [t 'lam P']);
t_is(r.branch(:, MU_SF) + r.branch(:, MU_ST), [0; 33; 0], 7, [t 'mu flow']);
t_is(r.reserves.R, [44; 100; 6; 0; 0], 7, [t 'R']);
t_is(r.reserves.prc, [5; 5; 5; 0; 0], 7, [t 'reserve prc']);
t_is(r.reserves.mu.Pmax + r.gen(:, MU_PMAX), [4; 0; 0; 0; 40], 7, [t 'reserve muPmax']);

%% most
t = sprintf('%s : DC OPF (w/reserves) : most   ', solvers{s});
mpc = mpc0;
% mpc.gen(4, PMIN) = -350.8;
mpopt = mpoption(mpopt, 'most.dc_model', 1);
mdi = loadmd(mpc, [], xgd);
mdi.FixedReserves = mpc.reserves;
mdo = most(mdi, mpopt);
rr = mdo.flow(1,1,1).mpc;
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -485656, 5, [t 'f']);
t_is(rr.gen(:, PG), [156; 65; 178; -499; 100], 6, [t 'Pg']);
t_is(rr.gen(:, GEN_STATUS), [1; 1; 1; 1; 1], 7, [t 'u']);
% rr.gen(:, GEN_STATUS)
t_is(rr.bus(:, LAM_P), [29; 40; 51], 6, [t 'lam P']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 33; 0], 7, [t 'mu flow']);
t_is(rr.reserves.R, [44; 100; 6; 0; 0], 6, [t 'R']);
t_is(rr.reserves.prc, [5; 5; 5; 0; 0], 7, [t 'reserve prc']);
t_is(rr.reserves.mu.Pmax + rr.gen(:, MU_PMAX), [4; 0; 0; 0; 40], 7, [t 'reserve muPmax']);

t = sprintf('%s : Secure DC OPF (w/cont,res,ramp) : c3sopf ', solvers{s});
mpc = mpc0;
mpc.gen(4, PMIN) = -350;
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
    1e-8    200     2e-8    0       1e-9    1e-9;
];
% xgd = loadxgendata(xgd_table, mpc);
xgd.PositiveActiveReserveQuantity(2) = 100;
xgd.PositiveActiveReserveQuantity(4) = 800;
xgd.NegativeActiveReserveQuantity(:) = 0;
xgd.PositiveActiveReservePrice([1;3]) = [5; 1.5];
xgd.NegativeActiveReservePrice([1;3]) = [10; 3];

mpc1 = mpc;
mpc1.gen(2, GEN_STATUS) = 0;
r = c3sopf(mpc1, xgd_table.data, contab, mpopt);
% r
% r.opf_results
% r.opf_results.f
t_ok(r.success, [t 'success']);
t_is(r.opf_results.f, -340350, 4, [t 'f']);
rr = r.base;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [190; 0; 60; -350; 100], 6, [t 'Pg base']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [25; 25; 25], 7, [t 'lam P base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow base']);
rr = r.cont(1);
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_ok(isempty(rr.gen), [t 'Pg 1']);
t_ok(isempty(rr.bus), [t 'lam P 1']);
t_ok(isempty(rr.branch), [t 'mu flow 1']);
rr = r.cont(2);
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [190; 0; 60; -300; 50], 5, [t 'Pg 2']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [0; 0; 50], 7, [t 'lam P 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 50], 7, [t 'mu flow 2']);

t = sprintf('%s : Secure DC OPF (w/cont,res,ramp) : most ', solvers{s});
mdi = loadmd(mpc, [], xgd, [], contab);
mdo = most(mdi, mpopt);
% mdo
% mdo.QP
% mdo.QP.f
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -340350, 5, [t 'f']);
rr = mdo.flow(1,1,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [190; 0; 60; -350; 100], 6, [t 'Pg base']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
% t_is(rr.bus(:, LAM_P), [20; 20; 20], 7, [t 'lam P base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow base']);
rr1 = rr;
rr = mdo.flow(1,1,2).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [190; 0; 60; -350; 100], 5, [t 'Pg 1']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(rr1.bus(:, LAM_P) + rr.bus(:, LAM_P), [25; 25; 25], 7, [t 'lam P base + lam P 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 7, [t 'mu flow 1']);
rr = mdo.flow(1,1,3).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [190; 0; 60; -300; 50], 5, [t 'Pg 2']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [0; 0; 50], 7, [t 'lam P 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 50], 7, [t 'mu flow 2']);

t = sprintf('%s : Stochastic DC OPF (w/wind,res) : c3sopf ', solvers{s});
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
mpc1 = mpc;
mpc1.gen(iwind, PMAX) = mpc1.gen(iwind, PMAX) * profiles.values(2);
mpc1.gen(3, GEN_STATUS) = 0;
wind = loadgenericdata('ex_wind', 'struct', 'gen', 'wind');
offers = [xgd_table.data; wind.xgd_table.data(:, [3:8])];
% mpopt.verbose = 2;
% mpopt.out.all = -1;
r = c3sopf(mpc1, offers, contab0, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.opf_results.f, -341928.605134, 5, [t 'f']);
rr = r.cont(1);
% rr.gen(:, PG)
% rr.bus(:, LAM_P)
% rr.branch(:, MU_SF) + rr.branch(:, MU_ST)
t_is(rr.gen(:, PG), [200; 150; 0; -350; 0], 7, [t 'Pg 1']);
t_is(rr.bus(:, LAM_P), [4.7596576; 4.7596576; 4.7596576], 7, [t 'lam P 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 7, [t 'mu flow 1']);
rr = r.base;
t_is(rr.gen(:, PG), [200; 100; 0; -350; 50], 5, [t 'Pg 2']);
t_is(rr.bus(:, LAM_P), [20.4806848; 20.4806848; 20.4806848], 6, [t 'lam P 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 2']);
rr = r.cont(2);
t_is(rr.gen(:, PG), [200; 65; 0; -350; 85], 5, [t 'Pg 3']);
t_is(rr.bus(:, LAM_P), [0; 0; 0], 6, [t 'lam P 3']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 7, [t 'mu flow 3']);

t = sprintf('%s : Stochastic DC OPF (w/wind,res) : most ', solvers{s});
mdi = loadmd(mpc, transmat, xgd, [], [], profiles);
mdo = most(mdi, mpopt);
% mdo
% mdo.QP
% mdo.QP.f
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -341928.605134, 6, [t 'f']);
rr = mdo.flow(1,1,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 150; 0; -350; 0], 7, [t 'Pg 1']);
t_is(rr.gen(:, GEN_STATUS), [1; 1; 0; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [4.7596576; 4.7596576; 4.7596576], 7, [t 'lam P 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 7, [t 'mu flow 1']);
rr = mdo.flow(1,2,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 100; 0; -350; 50], 7, [t 'Pg 2']);
t_is(rr.gen(:, GEN_STATUS), [1; 1; 0; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [20.4806848; 20.4806848; 20.4806848], 6, [t 'lam P 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 2']);
rr = mdo.flow(1,3,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 65; 0; -350; 85], 5, [t 'Pg 3']);
t_is(rr.gen(:, GEN_STATUS), [1; 1; 0; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [0; 0; 0], 6, [t 'lam P 3']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 7, [t 'mu flow 3']);
% keyboard;

t = sprintf('%s : Secure Stochastic DC OPF (w/wind,cont,res,ramp) : most ', solvers{s});
mdi = loadmd(mpc, transmat, xgd, [], contab, profiles);
mdo = most(mdi, mpopt);
% mdo
% mdo.QP
% mdo.QP.f
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -338372.01858, 4, [t 'f']);

rr = mdo.flow(1,1,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 0; 150; -350; 0], 7, [t 'Pg 1 base']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
% t_is(rr.bus(:, LAM_P), [6.2596576; 6.2596576; 6.2596576], 6, [t 'lam P 1 base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 2, [t 'mu flow 1 base']);

rr1 = rr;
rr = mdo.flow(1,1,2).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 0; 150; -350; 0], 4, [t 'Pg 1 1']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
% t_is(rr.bus(:, LAM_P), [1.2692420; 1.2692420; 1.2692420], 7, [t 'lam P 1 1']);
t_is(rr1.bus(:, LAM_P) + rr.bus(:, LAM_P), [7.5288996; 7.5288996; 7.5288996], 7, [t 'lam P 1 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 1 1']);

rr = mdo.flow(1,1,3).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 0; 100; -300; 0], 4, [t 'Pg 1 2']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [0.3173105; 0.3173105; 7.9327627], 6, [t 'lam P 1 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 7.6154522], 6, [t 'mu flow 1 2']);

rr = mdo.flow(1,2,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 0; 100; -350; 50], 6, [t 'Pg 2 base']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [20.4806848; 20.4806848; 20.4806848], 6, [t 'lam P 2 base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 2 base']);

rr = mdo.flow(1,2,2).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 0; 100; -350; 50], 6, [t 'Pg 2 1']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [5.4615159; 5.4615159; 5.4615159], 6, [t 'lam P 2 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 2 1']);

rr = mdo.flow(1,2,3).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 0; 60; -300; 40], 6, [t 'Pg 2 2']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [0; 0; 34.1344746], 5, [t 'lam P 2 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 34.1344746], 6, [t 'mu flow 2 2']);

rr = mdo.flow(1,3,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 0; 60; -350; 90], 6, [t 'Pg 3 base']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [0; 0; 0], 6, [t 'lam P 3 base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 5, [t 'mu flow 3 base']);

rr = mdo.flow(1,3,2).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 0; 60; -350; 90], 6, [t 'Pg 3 1']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [0; 0; 0], 6, [t 'lam P 3 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 3 1']);

rr = mdo.flow(1,3,3).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 0; 60; -300; 40], 4, [t 'Pg 3 2']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [0; 0; 7.9327627], 5, [t 'lam P 3 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 7.9327627], 6, [t 'mu flow 3 2']);
% keyboard;
end

end

t_end;
