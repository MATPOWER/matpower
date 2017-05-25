function t_most_spuc(quiet, create_plots, create_pdfs, savedir)
%T_MOST_SPUC  Tests of single-period unit commitment optimizations
%
%   T_MOST_SPUC(QUIET, CREATE_PLOTS, CREATE_PDFS, SAVEDIR)
%   Can generate summary plots and save them as PDFs in a directory of
%   your choice.
%   E.g. t_most_spuc(0, 1, 1, '~/Downloads/spuc_plots')

%   MOST
%   Copyright (c) 2015-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

if nargin < 4
    savedir = '.';              %% save in current working directory by default
    if nargin < 3
        create_pdfs = 0;        %% do NOT save plots to PDF files
        if nargin < 2
            create_plots = 0;   %% do NOT create summary plots of results
            if nargin < 1
                quiet = 0;      %% verbose by default
            end
        end
    end
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
    sc = mosek_symbcon;
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', sc.MSK_OPTIMIZER_FREE);            %% default
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', sc.MSK_OPTIMIZER_INTPNT);          %% interior point
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', sc.MSK_OPTIMIZER_PRIMAL_SIMPLEX);  %% primal simplex
    mpopt = mpoption(mpopt, 'mosek.lp_alg', sc.MSK_OPTIMIZER_DUAL_SIMPLEX);     %% dual simplex
    %mpopt = mpoption(mpopt, 'mosek.lp_alg', sc.MSK_OPTIMIZER_FREE_SIMPLEX);    %% automatic simplex
    %mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_MIO_TOL_X', 0);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_IPAR_MIO_NODE_OPTIMIZER', sc.MSK_OPTIMIZER_DUAL_SIMPLEX);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_IPAR_MIO_ROOT_OPTIMIZER', sc.MSK_OPTIMIZER_DUAL_SIMPLEX);
    mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_MIO_TOL_ABS_RELAX_INT', 1e-9);
    %mpopt = mpoption(mpopt, 'mosek.opts.MSK_DPAR_MIO_TOL_REL_RELAX_INT', 0);
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
    %% next line is to work around a bug in intlinprog
    % (Technical Support Case #01841662)
    mpopt = mpoption(mpopt, 'intlinprog.LPPreprocess', 'none');
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

xgd = loadxgendata('ex_xgd_uc', mpc);
[iwind, mpc, xgd] = addwind('ex_wind_uc', mpc, xgd);
mpc.reserves.zones = [mpc.reserves.zones 0];
mpc = scale_load(499, mpc, [], struct('scale', 'QUANTITY'));
mpc0 = mpc;
xgd0 = xgd;

%% data structures for results for plotting
if create_plots
    j = 1;
    Pg   = NaN(5, 7);
    Rp   = NaN(5, 7);
    Rm   = NaN(5, 7);
    lamP = NaN(3, 7);
    muF  = zeros(3, 7);
end

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
if create_plots
    Pg(:, j) = mdo.results.ExpectedDispatch;
    Rp(:, j) = 0;
    Rm(:, j) = 0;
    lamP(:, j) = rr.bus(:, LAM_P);
    j = j + 1;
end


%%-----  DC OPF  -----
if verbose
    fprintf('\n\n');
    fprintf('--------------------------------------------\n');
    fprintf('-----              DC OPF              -----\n');
    fprintf('--------------------------------------------\n');
end
t = sprintf('%s : DC OPF : runopf ', solvers{s});
mpc = mpc0;
% mpc = scale_load(500, mpc, [], struct('scale', 'QUANTITY'));
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
% mpc = scale_load(500, mpc, [], struct('scale', 'QUANTITY'));
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
if create_plots
    Pg(:, j) = mdo.results.ExpectedDispatch;
    Rp(:, j) = 0;
    Rm(:, j) = 0;
    lamP(:, j) = rr.bus(:, LAM_P);
    muF(:, j)  = rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
    j = j + 1;
end


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
% mpc = scale_load(350.8, mpc, [], struct('scale', 'QUANTITY'));
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
if create_plots
    Pg(:, j) = mdo.results.ExpectedDispatch;
    Rp(:, j) = rr.reserves.R;
    Rm(:, j) = 0;
    lamP(:, j) = rr.bus(:, LAM_P);
    j = j + 1;
end


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
% mpc = scale_load(350.8, mpc, [], struct('scale', 'QUANTITY'));
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
if create_plots
    Pg(:, j) = mdo.results.ExpectedDispatch;
    Rp(:, j) = rr.reserves.R;
    Rm(:, j) = 0;
    lamP(:, j) = rr.bus(:, LAM_P);
    muF(:, j)  = rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
    j = j + 1;
end

t = sprintf('%s : Secure DC OPF (w/cont,res,ramp) : c3sopf ', solvers{s});
mpc = mpc0;
mpc = scale_load(350, mpc, [], struct('scale', 'QUANTITY'));
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
t_is(r.opf_results.f, -340850, 4, [t 'f']);
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
t_is(rr.bus(:, LAM_P), [0; 0; 40], 7, [t 'lam P 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 40], 7, [t 'mu flow 2']);

t = sprintf('%s : Secure DC OPF (w/cont,res,ramp) : most ', solvers{s});
mdi = loadmd(mpc, [], xgd, [], contab);
mdo = most(mdi, mpopt);
% mdo
% mdo.QP
% mdo.QP.f
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -340850, 5, [t 'f']);
rr = mdo.flow(1,1,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [190; 0; 60; -350; 100], 6, [t 'Pg base']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
% t_is(rr.bus(:, LAM_P), [20; 20; 20], 7, [t 'lam P base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow base']);
if create_plots
    lamP(:, j) = rr.bus(:, LAM_P);
    muF(:, j)  = rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

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
if create_plots
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,1,3).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [190; 0; 60; -300; 50], 5, [t 'Pg 2']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [0; 0; 40], 7, [t 'lam P 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 40], 7, [t 'mu flow 2']);
if create_plots
    Pg(:, j) = mdo.results.ExpectedDispatch;
    Rp(:, j) = mdo.results.Rpp + mdo.results.Pc - mdo.results.ExpectedDispatch;
    Rm(:, j) = mdo.results.Rpm - mdo.results.Pc + mdo.results.ExpectedDispatch;
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
    j = j + 1;
end

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
wind = loadgenericdata('ex_wind_uc', 'struct', 'gen', 'wind');
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
if create_plots
    lamP(:, j) = rr.bus(:, LAM_P);
    muF(:, j)  = rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,2,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 100; 0; -350; 50], 7, [t 'Pg 2']);
t_is(rr.gen(:, GEN_STATUS), [1; 1; 0; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [20.4806848; 20.4806848; 20.4806848], 6, [t 'lam P 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 2']);
if create_plots
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,3,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 65; 0; -350; 85], 5, [t 'Pg 3']);
t_is(rr.gen(:, GEN_STATUS), [1; 1; 0; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [0; 0; 0], 6, [t 'lam P 3']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 7, [t 'mu flow 3']);
if create_plots
    Pg(:, j) = mdo.results.ExpectedDispatch;
    Rp(:, j) = mdo.results.Rpp + mdo.results.Pc - mdo.results.ExpectedDispatch;
    Rm(:, j) = mdo.results.Rpm - mdo.results.Pc + mdo.results.ExpectedDispatch;
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
    j = j + 1;
end
% keyboard;

t = sprintf('%s : Secure Stochastic DC OPF (w/wind,cont,res,ramp) : most ', solvers{s});
mdi = loadmd(mpc, transmat, xgd, [], contab, profiles);
mdo = most(mdi, mpopt);
% mdo
% mdo.QP
% mdo.QP.f
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -338857.9224447, 4, [t 'f']);
% fprintf('%.7f\n', mdo.QP.f);

rr = mdo.flow(1,1,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 0; 150; -350; 0], 7, [t 'Pg 1 base']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
% t_is(rr.bus(:, LAM_P), [7.2115891; 7.2115891; 7.2115891], 6, [t 'lam P 1 base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 2, [t 'mu flow 1 base']);
if create_plots
    lamP(:, j) = rr.bus(:, LAM_P);
    muF(:, j)  = rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr1 = rr;
rr = mdo.flow(1,1,2).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 0; 150; -350; 0], 4, [t 'Pg 1 1']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
% t_is(rr.bus(:, LAM_P), [0.3807726; 0.3807726; 0.3807726], 7, [t 'lam P 1 1']);
t_is(rr1.bus(:, LAM_P) + rr.bus(:, LAM_P), [7.5923617; 7.5923617; 7.5923617], 7, [t 'lam P 1 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 1 1']);
if create_plots
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,1,3).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 0; 100; -300; 0], 4, [t 'Pg 1 2']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [0.2538484; 0.2538484; 6.3462102], 6, [t 'lam P 1 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 6.0923618], 6, [t 'mu flow 1 2']);
if create_plots
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,2,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 0; 100; -350; 50], 6, [t 'Pg 2 base']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [24.5768217; 24.5768217; 24.5768217], 6, [t 'lam P 2 base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 2 base']);
if create_plots
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,2,2).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 0; 100; -350; 50], 6, [t 'Pg 2 1']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [1.6384548; 1.6384548; 1.6384548], 6, [t 'lam P 2 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 2 1']);
if create_plots
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,2,3).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 0; 60; -300; 40], 6, [t 'Pg 2 2']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [0; 0; 27.3075797], 5, [t 'lam P 2 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 27.3075797], 6, [t 'mu flow 2 2']);
if create_plots
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,3,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 0; 60; -350; 90], 6, [t 'Pg 3 base']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [0; 0; 0], 6, [t 'lam P 3 base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 5, [t 'mu flow 3 base']);
if create_plots
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,3,2).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 0; 60; -350; 90], 6, [t 'Pg 3 1']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [0; 0; 0], 6, [t 'lam P 3 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 3 1']);
if create_plots
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,3,3).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%d; %d; %d; %d; %d]\n', rr.gen(:, GEN_STATUS));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [200; 0; 60; -300; 40], 4, [t 'Pg 3 2']);
t_is(rr.gen(:, GEN_STATUS), [1; 0; 1; 1; 1], 7, [t 'u']);
t_is(rr.bus(:, LAM_P), [0; 0; 6.3462102], 5, [t 'lam P 3 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 6.3462102], 6, [t 'mu flow 3 2']);
% keyboard;
end

if create_plots
    create_plots = 0;   %% don't do them again

    Pg(:, j) = mdo.results.ExpectedDispatch;
    Rp(:, j) = mdo.results.Rpp + mdo.results.Pc - mdo.results.ExpectedDispatch;
    Rm(:, j) = mdo.results.Rpm - mdo.results.Pc + mdo.results.ExpectedDispatch;
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
%    R(4:5, :) = NaN;
    R = Rp + Rm;

    labels = {'Economic Dispatch'; 'DC OPF'; 'Economic Dispatch w/reserves'; 'DC OPF w/reserves'; 'secure DC OPF'; 'stochastic DC OPF'; 'secure stochastic DC OPF'};

    figure(1);

    subplot(4, 1, 1);
    bar(abs(Pg([1:3 5],:)'),'stacked');
    title('Generation');
    ylabel('MW');
    h = gca;
    h.XTickLabel = labels;
    h.XTickLabelRotation = 20;

    subplot(4, 1, 2);
    bar(abs(R([1:3 5],:)'),'stacked');
    title('Reserves');
    ylabel('MW');
    h = gca;
    h.XTickLabel = labels;
    h.XTickLabelRotation = 20;

    subplot(4, 1, 3);
    plot(lamP');
    title('Nodal Prices');
    ylabel('$/MW');
    h = gca;
    h.XTickLabel = {'', labels{:}, ''}';
    h.XTickLabelRotation = 20;
    h = [0 8 0 100];
    axis(h);

    subplot(4, 1, 4);
    plot(muF');
    title('Flow Constraint Shadow Prices');
    ylabel('$/MW');
    h = gca;
    h.XTickLabel = {'', labels{:}, ''}';
    h.XTickLabelRotation = 20;
    h = [0 8 0 60];
    axis(h);

    if create_pdfs
        fname = 'single-period-uc';
        h = gcf;
        set(h, 'PaperSize', [11 8.5]);
        set(h, 'PaperPosition', [0.25 0.25 10.5 8]);
        print('-dpdf', fullfile(savedir, fname));
    end

    for j = 1:7;
        figure(j+1);
        if create_pdfs
            fname = sprintf('single-period-uc-%d', j);
        else
            fname = '';
        end
        plot_case(labels{j}, Pg([1:3 5], j), Rp([1:3 5], j), Rm([1:3 5], j), lamP(:, j), muF(:, j), 250, 100, savedir, fname);
    end
end
end

t_end;


function h = plot_case(label, Pg, Rp, Rm, lamP, muF, maxq, maxp, mypath, fname)

subplot(1, 3, 1);
h = bar([Pg-Rm Rm Rp], 'stacked');
set(h(2), 'FaceColor', [0 0.35 0.33]);
ah1 = gca;
title('Generation & Reserves');
ylabel('MW');
xlabel('Gen');
set(gca, 'FontSize', 14);

if nargin < 6
    maxq = ah1.YLim(2);
end
ah1.YLim(2) = maxq;
ah1.YLim(1) = 0;

subplot(1, 3, 2);
bar(lamP);
ah3 = gca;
title('Nodal Prices');
ylabel('$/MW');
xlabel('Bus');
set(gca, 'FontSize', 14);

subplot(1, 3, 3);
bar(muF);
ah4 = gca;
title('Flow Constraint Shadow Prices');
ylabel('$/MW');
xlabel('Line');
set(gca, 'FontSize', 14);

if nargin < 7
    maxp = max(ah3.YLim(2), ah4.YLim(2));
end
ah3.YLim(1) = 0;
ah4.YLim(1) = 0;
ah3.YLim(2) = maxp;
ah4.YLim(2) = maxp;

[ax,h] = suplabel(label, 't');
set(h, 'FontSize', 18)

if nargin > 7 && ~isempty(fname)
    h = gcf;
    set(h, 'PaperSize', [11 8.5]);
    set(h, 'PaperPosition', [0.25 0.25 10.5 8]);
    print('-dpdf', fullfile(mypath, fname));
end
