function t_most_sp(quiet, create_plots, create_pdfs, savedir)
%T_MOST_SP  Tests of single-period continuous optimizations
%
%   T_MOST_SP(QUIET, CREATE_PLOTS, CREATE_PDFS, SAVEDIR)
%   Can generate summary plots and save them as PDFs in a directory of
%   your choice.
%   E.g. t_most_sp(0, 1, 1, '~/Downloads/sp_plots')

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

t_begin(156, quiet);

if quiet
    verbose = 0;
else
    verbose = 0;
end
% verbose = 2;

casefile = 'ex_case3a';
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

%% turn off warnings
if have_fcn('octave')
    s = warning('query', 'Octave:nearly-singular-matrix');
    warning('off', 'Octave:nearly-singular-matrix');
end

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
xgd_table = loadgenericdata('ex_xgd', 'struct');

%%-----  wind  -----
wind = loadgenericdata('ex_wind', 'struct');

%%-----  contingencies  -----
contab = loadgenericdata('ex_contab', 'array');
pp = [1-sum(contab(:,2)); contab(1,2); contab(2,2)];

xgd = loadxgendata(xgd_table, mpc);

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
mdi = loadmd(mpc);
mdo = most(mdi, mpopt);
rr = mdo.flow(1,1,1).mpc;
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -443250, 5, [t 'f']);
t_is(rr.gen(:, PG), [150; 150; 150; -450], 7, [t 'Pg']);
t_is(rr.bus(:, LAM_P), [30; 30; 30], 7, [t 'lam P']);
if create_plots
    Pg(1:4, j) = mdo.results.ExpectedDispatch;
    Rp(1:4, j) = 0;
    Rm(1:4, j) = 0;
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
mdi = loadmd(mpc);
mdo = most(mdi, mpopt);
rr = mdo.flow(1,1,1).mpc;
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -443115, 5, [t 'f']);
t_is(rr.gen(:, PG), [135; 135; 180; -450], 6, [t 'Pg']);
t_is(rr.bus(:, LAM_P), [27; 36; 45], 7, [t 'lam P']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 27; 0], 6, [t 'mu flow']);
if create_plots
    Pg(1:4, j) = mdo.results.ExpectedDispatch;
    Rp(1:4, j) = 0;
    Rm(1:4, j) = 0;
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
mdi = loadmd(mpc);
mdi.FixedReserves = mpc.reserves;
mdo = most(mdi, mpopt);
rr = mdo.flow(1,1,1).mpc;
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -442820, 5, [t 'f']);
t_is(rr.gen(:, PG), [140; 150; 160; -450], 7, [t 'Pg']);
t_is(rr.bus(:, LAM_P), [32; 32; 32], 7, [t 'lam P']);
t_is(rr.reserves.R, [60; 50; 40; 0], 7, [t 'R']);
t_is(rr.reserves.prc, [5; 5; 5; 0], 7, [t 'reserve prc']);
t_is(rr.reserves.mu.Pmax + rr.gen(:, MU_PMAX), [4; 2; 0; 0], 7, [t 'reserve muPmax']);
if create_plots
    Pg(1:4, j) = mdo.results.ExpectedDispatch;
    Rp(1:4, j) = rr.reserves.R;
    Rm(1:4, j) = 0;
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
mdi = loadmd(mpc);
mdi.FixedReserves = mpc.reserves;
mdo = most(mdi, mpopt);
rr = mdo.flow(1,1,1).mpc;
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -442760, 5, [t 'f']);
t_is(rr.gen(:, PG), [130; 140; 180; -450], 6, [t 'Pg']);
t_is(rr.bus(:, LAM_P), [30; 36; 42], 6, [t 'lam P']);
t_is(rr.reserves.R, [70; 60; 20; 0], 6, [t 'R']);
t_is(rr.reserves.prc, [5; 5; 5; 0], 7, [t 'reserve prc']);
t_is(rr.reserves.mu.Pmax + rr.gen(:, MU_PMAX), [4; 2; 0; 0], 7, [t 'reserve muPmax']);
if create_plots
    Pg(1:4, j) = mdo.results.ExpectedDispatch;
    Rp(1:4, j) = rr.reserves.R;
    Rm(1:4, j) = 0;
    lamP(:, j) = rr.bus(:, LAM_P);
    muF(:, j)  = rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
    j = j + 1;
end


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
t_is(mdo.QP.f, sum(pp.*ff), 3, [t 'f']);
rr = mdo.flow(1,1,1).mpc;
t_is(rr.gen(:, PG), [135; 135; 180; -450], 5, [t 'Pg base']);
t_is(rr.bus(:, LAM_P), [27; 36; 45]*pp(1), 5, [t 'lam P base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 27; 0]*pp(1), 5, [t 'mu flow base']);
rr = mdo.flow(1,1,2).mpc;
t_is(rr.gen(:, PG), [200; 0; 250; -450], 4, [t 'Pg 1']);
t_is(rr.bus(:, LAM_P), [50; 50; 50]*pp(2), 5, [t 'lam P 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0]*pp(2), 5, [t 'mu flow 1']);
rr = mdo.flow(1,1,3).mpc;
t_is(rr.gen(:, PG), [100; 100; 100; -300], 4, [t 'Pg 2']);
t_is(rr.bus(:, LAM_P), [20; 20; 1000]*pp(3), 6, [t 'lam P 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 980]*pp(3), 6, [t 'mu flow 2']);

% mpopt.verbose = 2;
% mpopt.out.all = -1;

t = 'Secure DC OPF (w/cont,res,ramp) : c3sopf ';
xgd_table = loadgenericdata('ex_xgd_res', 'struct');
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
t_is(r.opf_results.f, -436686.1245, 5, [t 'f']);
rr = r.base;
% rr.gen(:, PG)
% rr.bus(:, LAM_P)
% rr.branch(:, MU_SF) + rr.branch(:, MU_ST)
t_is(rr.gen(:, PG), [142.15; 127.85; 180; -450], 7, [t 'Pg base']);
t_is(rr.bus(:, LAM_P), [23.6958; 32.4; 41.1042], 7, [t 'lam P base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 26.1126; 0], 7, [t 'mu flow base']);
rr = r.cont(1);
t_is(rr.gen(:, PG), [142.15; 0; 307.85; -450], 6, [t 'Pg 1']);
t_is(rr.bus(:, LAM_P), [5.1942; 5.1942; 5.1942], 7, [t 'lam P 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 7, [t 'mu flow 1']);
rr = r.cont(2);
t_is(rr.gen(:, PG), [142.15; 27.85; 130; -300], 5, [t 'Pg 2']);
t_is(rr.bus(:, LAM_P), [-0.46; -0.46; 40], 7, [t 'lam P 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 40.46], 7, [t 'mu flow 2']);

t = 'Secure DC OPF (w/cont,res,ramp) : most ';
mdi = loadmd(mpc, [], xgd, [], contab);
mdo = most(mdi, mpopt);
% mdo
% mdo.QP
% mdo.QP.f
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -436686.1245, 2, [t 'f']);
rr = mdo.flow(1,1,1).mpc;
t_is(rr.gen(:, PG), [142.15; 127.85; 180; -450], 4, [t 'Pg base']);
t_is(rr.bus(:, LAM_P), [23.6958; 32.4; 41.1042], 5, [t 'lam P base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 26.1126; 0], 5, [t 'mu flow base']);
if create_plots
    lamP(:, j) = rr.bus(:, LAM_P);
    muF(:, j)  = rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,1,2).mpc;
t_is(rr.gen(:, PG), [142.15; 0; 307.85; -450], 4, [t 'Pg 1']);
t_is(rr.bus(:, LAM_P), [5.1942; 5.1942; 5.1942], 6, [t 'lam P 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 1']);
if create_plots
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,1,3).mpc;
t_is(rr.gen(:, PG), [142.15; 27.85; 130; -300], 3, [t 'Pg 2']);
t_is(rr.bus(:, LAM_P), [-0.46; -0.46; 40], 6, [t 'lam P 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 40.46], 6, [t 'mu flow 2']);
if create_plots
    Pg(1:4, j) = mdo.results.ExpectedDispatch;
    Rp(1:4, j) = mdo.results.Rpp + mdo.results.Pc - mdo.results.ExpectedDispatch;
    Rm(1:4, j) = mdo.results.Rpm - mdo.results.Pc + mdo.results.ExpectedDispatch;
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
    j = j + 1;
end

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
if create_plots
    lamP(:, j) = rr.bus(:, LAM_P);
    muF(:, j)  = rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,2,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [128.7823267; 135.6088362; 135.6088370; -450; 50], 5, [t 'Pg 2']);
t_is(rr.bus(:, LAM_P), [18.5157455; 18.5157455; 18.5157455], 6, [t 'lam P 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 2']);
if create_plots
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,3,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [128.7823267; 86.9726845; 134.2449888; -450; 100], 5, [t 'Pg 3']);
t_is(rr.bus(:, LAM_P), [2.7597346; 2.7597346; 2.7597346], 6, [t 'lam P 3']);
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


t = 'Secure Stochastic DC OPF (w/wind,cont,res,ramp) : most ';
mdi = loadmd(mpc, transmat, xgd, [], contab, profiles);
mdo = most(mdi, mpopt);
% mdo
% mdo.QP
% mdo.QP.f
t_ok(mdo.QP.exitflag > 0, [t 'success']);
t_is(mdo.QP.f, -438182.429, 5, [t 'f']);

rr = mdo.flow(1,1,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [137.6930222; 132.3069578; 180; -450; 0], 4, [t 'Pg 1 base']);
t_is(rr.bus(:, LAM_P), [3.9958603; 5.1404304; 6.2850005], 3, [t 'lam P 1 base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 3.4337103; 0], 2, [t 'mu flow 1 base']);
if create_plots
    lamP(:, j) = rr.bus(:, LAM_P);
    muF(:, j)  = rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,1,2).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [137.6930928; 0; 312.3069041; -450; 0], 5, [t 'Pg 1 1']);
t_is(rr.bus(:, LAM_P), [2.0945890; 2.0945890; 2.0945892], 6, [t 'lam P 1 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 1 1']);
if create_plots
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,1,3).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [137.6930139; 45.2766586; 117.0303189; -300; 0], 4, [t 'Pg 1 2']);
t_is(rr.bus(:, LAM_P), [0.0574664; 0.0574664; 6.3462102], 6, [t 'lam P 1 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 6.2887437], 6, [t 'mu flow 1 2']);
if create_plots
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,2,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [137.6929448; 131.1533869; 131.1535795; -450; 50], 4, [t 'Pg 2 base']);
t_is(rr.bus(:, LAM_P), [16.1166803; 16.1166890; 16.1166978], 6, [t 'lam P 2 base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 4, [t 'mu flow 2 base']);
if create_plots
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,2,2).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [137.6930711; 0; 262.3068485; -450; 50], 4, [t 'Pg 2 1']);
t_is(rr.bus(:, LAM_P), [2.1488905; 2.1488905; 2.1488907], 6, [t 'lam P 2 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 2 1']);
if create_plots
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,2,3).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [137.6929611; 32.3066174; 117.0299763; -300; 12.9704443], 6, [t 'Pg 2 2']);
t_is(rr.bus(:, LAM_P), [0; 0; 27.3075797], 6, [t 'lam P 2 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 27.3075797], 6, [t 'mu flow 2 2']);
if create_plots
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,3,1).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [137.6929341; 95.2771176; 117.0299561; -450; 100], 4, [t 'Pg 3 base']);
t_is(rr.bus(:, LAM_P), [2.7209142; 2.7209144; 2.7209147], 6, [t 'lam P 3 base']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 5, [t 'mu flow 3 base']);
if create_plots
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,3,2).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [137.6930070; 0.0000000; 212.3070512; -450; 100], 4, [t 'Pg 3 1']);
t_is(rr.bus(:, LAM_P), [0.4042034; 0.4042034; 0.4042036], 6, [t 'lam P 3 1']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 0], 6, [t 'mu flow 3 1']);
if create_plots
    lamP(:, j) = lamP(:, j) + rr.bus(:, LAM_P);
    muF(:, j)  = muF(:, j) + rr.branch(:, MU_SF) + rr.branch(:, MU_ST);
end

rr = mdo.flow(1,3,3).mpc;
% fprintf('[%.7f; %.7f; %.7f; %.7f; %.7f]\n', rr.gen(:, PG));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.bus(:, LAM_P));
% fprintf('[%.7f; %.7f; %.7f]\n', rr.branch(:, MU_SF) + rr.branch(:, MU_ST));
t_is(rr.gen(:, PG), [137.6930131; 32.3071796; 117.0302086; -300; 12.9695914], 5, [t 'Pg 3 2']);
t_is(rr.bus(:, LAM_P), [0; 0; 6.3462102], 6, [t 'lam P 3 2']);
t_is(rr.branch(:, MU_SF) + rr.branch(:, MU_ST), [0; 0; 6.3462102], 6, [t 'mu flow 3 2']);
% keyboard;

if create_plots
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
        fname = 'single-period-continuous';
        h = gcf;
        set(h, 'PaperSize', [11 8.5]);
        set(h, 'PaperPosition', [0.25 0.25 10.5 8]);
        print('-dpdf', fullfile(savedir, fname));
    end

    for j = 1:7;
        figure(j+1);
        if create_pdfs
            fname = sprintf('single-period-cont-%d', j);
        else
            fname = '';
        end
        plot_case(labels{j}, Pg([1:3 5], j), Rp([1:3 5], j), Rm([1:3 5], j), lamP(:, j), muF(:, j), 320, 100, savedir, fname);
    end
end

%% turn warnings back on
if have_fcn('octave')
    warning(s.state, 'Octave:nearly-singular-matrix');
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
