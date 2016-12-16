function t_most_3b_1_1_0(quiet)
%T_MOST_3B_1_1_0  Tests for MOST.

%   MOST
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

if nargin < 1
    quiet = 0;
end

n_tests = 15;

t_begin(n_tests, quiet);

casename = 't_case3_most';
fudging = struct( ...       %% paramters for fudging reserve contract for sopf2
    'fudge',    0.05, ...   %% initial value (MW)
    'step',     0.01, ...   %% if necessary, increase by this amount and retry (MW)
    'lim',      0.1);       %% upper limit (MW), give up if no convergence
                            %% with fudge equal to this limit

%% options
mpopt = mpoption('verbose', 0, 'out.all', 0);
mpopt = mpoption(mpopt, 'opf.violation', 5e-7, 'mips.comptol', 5e-8);
mpoptac = mpoption(mpopt, 'model', 'AC');
mpoptdc = mpoption(mpopt, 'model', 'DC');
mpopt = mpoption(mpopt, 'most.solver', 'DEFAULT');

%% turn off warnings
s7 = warning('query', 'MATLAB:nearlySingularMatrix');
s6 = warning('query', 'MATLAB:nearlySingularMatrixUMFPACK');
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:nearlySingularMatrixUMFPACK');

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
    1       400     2       400     0.01    0.01;
    3       300     4       300     0.01    0.01;
    0.001   450     0.002   450     0       0;
];

%% contingency table
% label probty  type        row column      chgtype newvalue
contab = [
%    1   0.01    CT_TBRCH    1   BR_STATUS   CT_REP  0;      %% line 1-2
%    2   0.01    CT_TGEN     1   GEN_STATUS  CT_REP  0;      %% gen 1 at bus 1
%   2   0.00049621965   CT_TGEN     1   GEN_STATUS  CT_REP  0;      %% gen 1 at bus 1
];
clist = [];
nc = length(clist);

%% load the case
mpc = loadcase(casename);
gbus = mpc.gen(:, GEN_BUS);

%%-----  get OPF results  -----
rdc = rundcopf(mpc, mpoptdc);
% rac = runopf(mpc, mpoptac);
% save t_most5_soln rdc rac -v6
% s = load('t_most5_soln');
s.rdc = rdc;
% s.rac = rac;

%%-----  set up data for DC run (most)  -----
ng = size(mpc.gen, 1);      %% number of gens
xgd = loadxgendata(xgd_table, mpc);
md = loadmd(mpc, [], xgd);

%%-----  do DC run (most)  -----
r = most(md, mpopt);

%%-----  test the results  -----
t = 'success1';
t_ok(s.rdc.success, t);
t = 'success2';
t_ok(r.QP.exitflag, t);

t = 'f';
t_is(r.results.f, s.rdc.f, 4, t);

t = 'Pg : base';
t_is(r.flow(1,1,1).mpc.gen(:, PG), s.rdc.gen(:, PG), 5, t);

t = 'gen : base';
t_is(r.flow(1,1,1).mpc.gen(:,1:MU_PMIN), s.rdc.gen(:,1:MU_PMIN), 3, t);

t = 'energy prices';
t_is(r.results.GenPrices, s.rdc.bus(gbus, LAM_P), 6, t);

t = 'Pc';
t_is(r.results.Pc, s.rdc.gen(:, PG), 4, t);

t = 'Gmin';
t_is(r.results.Pc - r.results.Rpm, s.rdc.gen(:, PG), 4, t);

t = 'Gmax';
t_is(r.results.Pc + r.results.Rpp, s.rdc.gen(:, PG), 4, t);

t = 'upward contingency reserve quantities';
t_is(r.results.Rpp, zeros(ng, 1), 4, t);

t = 'downward contingency reserve quantities';
t_is(r.results.Rpm, zeros(ng, 1), 4, t);

t = 'upward contingency reserve prices';
t_is(r.results.RppPrices, xgd.PositiveActiveReservePrice, 6, t);

t = 'downward contingency reserve prices';
t_is(r.results.RpmPrices, xgd.NegativeActiveReservePrice, 6, t);

t = 'Rpmax_pos';
vv = get_idx(r.om);
Rpmax_pos = (r.QP.lambda.upper(vv.i1.Rpp(1):vv.iN.Rpp(1)) - r.QP.lambda.lower(vv.i1.Rpp(1):vv.iN.Rpp(1))) / mpc.baseMVA;
t_is(Rpmax_pos, zeros(ng, 1), 6, t);

t = 'Rpmax_neg';
Rpmax_neg = (r.QP.lambda.upper(vv.i1.Rpm(1):vv.iN.Rpm(1)) - r.QP.lambda.lower(vv.i1.Rpm(1):vv.iN.Rpm(1))) / mpc.baseMVA;
t_is(Rpmax_neg, zeros(ng, 1), 6, t);



% g1 = s.rdc.base.gen(:, PG);
% g2 = r.flow(1,1,1).mpc.gen(:, PG);
% for k = 1:nc
%     g1 = [ g1 s.rdc.cont(k).gen(:, PG) ];
%     g2 = [ g2 r.flow(1,1,k+1).mpc.gen(:, PG) ];
% end
% [m,n] = size(g1);
% for j = 1:n
%     fprintf('\n');
%     for i = 1:m
%         fprintf('%9.2f  %9.2f\n', g1(i,j), g2(i,j));
%     end
% end

%%-----  do AC run (most)  -----
%mostac;




%% turn warnings back on
warning(s7.state, 'MATLAB:nearlySingularMatrix');
warning(s6.state, 'MATLAB:nearlySingularMatrixUMFPACK');

t_end;
