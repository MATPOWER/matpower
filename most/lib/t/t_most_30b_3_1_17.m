function t_most_30b_3_1_17(quiet)
%T_MOST_30B_3_1_17  Tests for MOST.

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

n_tests = 96;

t_begin(n_tests, quiet);

casename = 't_case30_most';
fudging = struct( ...       %% paramters for fudging reserve contract for sopf2
    'fudge',    0.05, ...   %% initial value (MW)
    'step',     0.01, ...   %% if necessary, increase by this amount and retry (MW)
    'lim',      0.1);       %% upper limit (MW), give up if no convergence
                            %% with fudge equal to this limit

%% options
algs.dc     = {'DEFAULT'};  %% opf.dc.solver sequence to try for c3sopf (DC run)
algs.ac     = {'DEFAULT'};  %% opf.ac.solver sequence to try for c3sopf (AC run)
mpopt = mpoption('verbose', 0, 'out.all', 0);
mpopt = mpoption(mpopt, 'opf.violation', 5e-7, 'mips.comptol', 5e-8);
mpopt = mpoption(mpopt, 'sopf.force_Pc_eq_P0', 0);  %% don't constrain contracted == base case dispatch
if have_fcn('gurobi')
    mpopt = mpoption(mpopt, 'gurobi.method', 1);    %% dual-simplex
end
if have_fcn('mosek')
    sc = mosek_symbcon;
    mpopt = mpoption(mpopt, 'mosek.lp_alg', sc.MSK_OPTIMIZER_DUAL_SIMPLEX);     %% dual simplex
end
if have_fcn('linprog')
    if have_fcn('linprog_ds')
        mpopt = mpoption(mpopt, 'linprog.Algorithm', 'dual-simplex');
    else
        mpopt = mpoption(mpopt, 'linprog.Algorithm', 'simplex');
    end
end
mpoptac = mpoption(mpopt, 'model', 'AC');
mpoptdc = mpoption(mpopt, 'model', 'DC');
mpopt = mpoption(mpopt, 'most.solver', algs.dc{1});

%% turn off warnings
s7 = warning('query', 'MATLAB:nearlySingularMatrix');
s6 = warning('query', 'MATLAB:nearlySingularMatrixUMFPACK');
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:nearlySingularMatrixUMFPACK');

%% define named indices into data matrices
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
    10.1    15      10.0    15      0.1     0.0;
    10.3    30      10.2    30      0.3     0.2;  
    10.5    20      10.4    20      0.5     0.4;  
    10.7    25      10.6    25      0.7     0.6;
    20.1    25      20.0    25      60.1    60.0; 
    20.3    15      20.2    15      15.1    15.0;
    20.5    30      20.4    30      60.3    60.2;
    20.7    15      20.6    15      15.3    15.2;
    30.1    15      30.0    15      60.3    60.4;
    30.3    30      30.2    30      30.1    30.0;
    30.5    25      30.4    25      60.7    60.6;
    30.7    30      30.6    30      30.3    30.2;
    0.001   50      0.002   50      0       0;
    0.001   50      0.002   50      0       0;
    0.001   50      0.002   50      0       0;
    0.001   50      0.002   50      0       0;
    0.001   50      0.002   50      0       0;
    0.001   50      0.002   50      0       0;
    0.001   50      0.002   50      0       0;
    0.001   50      0.002   50      0       0;
    0.001   50      0.002   50      0       0;
    0.001   50      0.002   50      0       0;
    0.001   50      0.002   50      0       0;
    0.001   50      0.002   50      0       0;
    0.001   50      0.002   50      0       0;
    0.001   50      0.002   50      0       0;
    0.001   50      0.002   50      0       0;
    0.001   50      0.002   50      0       0;
    0.001   50      0.002   50      0       0;
    0.001   50      0.002   50      0       0;
    0.001   50      0.002   50      0       0;
    0.001   50      0.002   50      0       0;
];

%% contingency table
% label probty  type        row column      chgtype newvalue
contab = [
    1   0.002   CT_TBRCH    1   BR_STATUS   CT_REP  0;      %% line 1-2
    2   0.002   CT_TBRCH    2   BR_STATUS   CT_REP  0;      %% line 1-3, all power from gen1 flows via gen2
    3   0.002   CT_TBRCH    3   BR_STATUS   CT_REP  0;      %% line 2-4, a path to loads @ buses 7 & 8
    4   0.002   CT_TBRCH    5   BR_STATUS   CT_REP  0;      %% line 2-5, a path to loads @ buses 7 & 8
    5   0.002   CT_TBRCH    6   BR_STATUS   CT_REP  0;      %% line 2-6, a path to loads @ buses 7 & 8
    6   0.002   CT_TBRCH    36  BR_STATUS   CT_REP  0;      %% line 28-27, tie line between areas 1 & 3
    7   0.002   CT_TBRCH    15  BR_STATUS   CT_REP  0;      %% line 4-12, tie line between areas 1 & 2
    8   0.002   CT_TBRCH    12  BR_STATUS   CT_REP  0;      %% line 6-10, tie line between areas 1 & 3
    9   0.002   CT_TBRCH    14  BR_STATUS   CT_REP  0;      %% line 9-10, tie line between areas 1 & 3
    10  0.002   CT_TGEN     1   GEN_STATUS  CT_REP  0;      %% gen 1 at bus 1             
    11  0.002   CT_TGEN     2   GEN_STATUS  CT_REP  0;      %% gen 2 at bus 2
    12  0.002   CT_TGEN     3   GEN_STATUS  CT_REP  0;      %% gen 3 at bus 22 
    13  0.002   CT_TGEN     4   GEN_STATUS  CT_REP  0;      %% gen 4 at bus 27
    14  0.002   CT_TGEN     5   GEN_STATUS  CT_REP  0;      %% gen 5 at bus 23
    15  0.002   CT_TGEN     6   GEN_STATUS  CT_REP  0;      %% gen 6 at bus 13
    20  0.010   CT_TGEN     0   PMIN        CT_REL  1.1;    %% 10% load increase
    20  0.010   CT_TGEN     0   QMIN        CT_REL  1.1;
    21  0.010   CT_TGEN     0   PMIN        CT_REL  0.9;    %% 10% load decrease
    21  0.010   CT_TGEN     0   QMIN        CT_REL  0.9;
];
clist = unique(contab(:, CT_LABEL));
nc = length(clist);

%% load the case
mpc = loadcase(casename);
gbus = mpc.gen(:, GEN_BUS);

%%-----  get c3sopf results  -----
rdc = c3sopf_retry(algs.dc, mpc, xgd_table.data, contab, mpoptdc);
% rac = c3sopf_retry(algs.ac, mpc, xgd_table.data, contab, mpoptac);
% save t_most2_soln rdc rac -v6
% s = load('t_most2_soln');
s.rdc = rdc;
% s.rac = rac;

%%-----  set up data for DC run (most)  -----
ng = size(mpc.gen, 1);      %% number of gens
nt = 3;
xgd = loadxgendata(xgd_table, mpc);
md = loadmd(mpc, nt, xgd, [], contab);

%%-----  do DC run (most)  -----
r = most(md, mpopt);

%%-----  test the results  -----
t = 'success1';
t_ok(s.rdc.opf_results.success, t);
t = 'success2';
t_ok(r.QP.exitflag, t);

t = 'f';
t_is(r.results.f/sum(r.StepProb), s.rdc.opf_results.f, 4, t);

for tt = 1:nt

t = sprintf('(t=%d) Pg : base', tt);
t_is(r.flow(tt,1,1).mpc.gen(:, PG), s.rdc.base.gen(:, PG), 5, t);
t = sprintf('(t=%d) Pg : cont ', tt);
for k = 1:nc
    t_is(r.flow(tt,1,k+1).mpc.gen(:, PG), s.rdc.cont(k).gen(:, PG), 5, sprintf('%s %d', t, k));
end

% t = sprintf('(t=%d) gen : base', tt);
% t_is(r.flow(tt,1,1).mpc.gen(:,1:MU_PMIN), s.rdc.base.gen(:,1:MU_PMIN), 3, t);
% t = sprintf('(t=%d) gen : cont ', tt);
% for k = 1:nc
%     t_is(r.flow(tt,1,k+1).mpc.gen(:,1:MU_PMIN), s.rdc.cont(k).gen(:,1:MU_PMIN), 3, sprintf('%s %d', t, k));
% end

t = sprintf('(t=%d) energy prices', tt);
t_is(r.results.GenPrices(:,tt)/r.StepProb(tt), s.rdc.energy.prc.sum_bus_lam_p(gbus), 6, t);

t = sprintf('(t=%d) Pc', tt);
t_is(r.results.Pc(:,tt), s.rdc.energy.Pc, 4, t);

t = sprintf('(t=%d) Gmin', tt);
t_is(r.results.Pc(:,tt) - r.results.Rpm(:,tt), s.rdc.energy.Gmin, 4, t);

t = sprintf('(t=%d) Gmax', tt);
t_is(r.results.Pc(:,tt) + r.results.Rpp(:,tt), s.rdc.energy.Gmax, 4, t);

t = sprintf('(t=%d) upward contingency reserve quantities', tt);
t_is(r.results.Rpp(:,tt), s.rdc.reserve.qty.Rp_pos, 4, t);

t = sprintf('(t=%d) downward contingency reserve quantities', tt);
t_is(r.results.Rpm(:,tt), s.rdc.reserve.qty.Rp_neg, 4, t);

t = sprintf('(t=%d) upward contingency reserve prices', tt);
t_is(r.results.RppPrices(:,tt)/r.StepProb(tt), s.rdc.reserve.prc.Rp_pos, 6, t);

t = sprintf('(t=%d) downward contingency reserve prices', tt);
t_is(r.results.RpmPrices(:,tt)/r.StepProb(tt), s.rdc.reserve.prc.Rp_neg, 6, t);

t = sprintf('(t=%d) contingency physical ramp price', tt);
[vv, ll] = get_idx(r.om);
Ramp_P_max = zeros(ng, nc);
sum_muPmax = zeros(ng, 1);
sum_muPmin = zeros(ng, 1);
for k = 1:nc+1
    ii = find(r.flow(tt,1,k).mpc.gen(:, GEN_STATUS) > 0);
    if k > 1
        Ramp_P_max(ii,k-1) = (r.QP.lambda.mu_u(ll.i1.rampcont(tt,1,k):ll.iN.rampcont(tt,1,k)) - r.QP.lambda.mu_l(ll.i1.rampcont(tt,1,k):ll.iN.rampcont(tt,1,k))) / mpc.baseMVA;
    end
    sum_muPmax(ii) = sum_muPmax(ii) + r.flow(tt,1,k).mpc.gen(ii, MU_PMAX);
    sum_muPmin(ii) = sum_muPmin(ii) + r.flow(tt,1,k).mpc.gen(ii, MU_PMIN);
end
t_is(sum(Ramp_P_max/r.StepProb(tt), 2), sum(s.rdc.energy.mu.Ramp_P_max, 2), 6, t);
% t_is(Ramp_P_max/r.StepProb(tt), s.rdc.energy.mu.Ramp_P_max, 6, t);

t = sprintf('(t=%d) sum_muPmax', tt);
t_is(sum_muPmax/r.StepProb(tt), s.rdc.energy.sum_muPmax, 1, t);

t = sprintf('(t=%d) sum_muPmin', tt);
t_is(sum_muPmin/r.StepProb(tt), s.rdc.energy.sum_muPmin, 0, t);

t = sprintf('(t=%d) Rpmax_pos', tt);
Rpmax_pos = (r.QP.lambda.upper(vv.i1.Rpp(1):vv.iN.Rpp(1)) - r.QP.lambda.lower(vv.i1.Rpp(1):vv.iN.Rpp(1))) / mpc.baseMVA;
t_is(Rpmax_pos, s.rdc.reserve.mu.Rpmax_pos, 6, t);

t = sprintf('(t=%d) Rpmax_neg', tt);
Rpmax_neg = (r.QP.lambda.upper(vv.i1.Rpm(1):vv.iN.Rpm(1)) - r.QP.lambda.lower(vv.i1.Rpm(1):vv.iN.Rpm(1))) / mpc.baseMVA;
t_is(Rpmax_neg, s.rdc.reserve.mu.Rpmax_neg, 6, t);

end

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
