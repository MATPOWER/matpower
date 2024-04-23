function t_makePTDF(quiet)
% t_makePTDF - Tests for makePTDF.

%   MATPOWER
%   Copyright (c) 2006-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

ntests = 67;
t_begin(ntests, quiet);

casefile = 't_case9_opf';
if quiet
    verbose = 0;
else
    verbose = 0;
end

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% load case
mpopt = mpoption('out.all', 0, 'verbose', verbose);
r = rundcopf(casefile, mpopt);
mpc = ext2int(r);
[baseMVA, bus, gen, branch] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);
% [baseMVA, bus, gen, gencost, branch, f, success, et] = ...
%     rundcopf(casefile, mpopt);
% [i2e, bus, gen, branch] = ext2int(bus, gen, branch);
nb  = size(bus, 1);
nbr = size(branch, 1);
ng  = size(gen, 1);

%% compute injections and flows
Cg = sparse(gen(:, GEN_BUS), (1:ng)', ones(ng, 1), nb, ng);
Pg = Cg * gen(:, PG);
Pd = bus(:, PD);
P  = Pg - Pd;
ig = find(P > 0);
il = find(P <= 0);
F  = branch(:, PF);

%% create corresponding slack distribution matrices
e1 = zeros(nb, 1);  e1(1) = 1;
e4 = zeros(nb, 1);  e4(4) = 1;
D1  = eye(nb, nb) - e1 * ones(1, nb);
D4  = eye(nb, nb) - e4 * ones(1, nb);
Deq = eye(nb, nb) - ones(nb, 1) / nb * ones(1, nb);
Dd  = eye(nb) - Pd/sum(Pd) * ones(1, nb);
Dg  = eye(nb) - Pg/sum(Pg) * ones(1, nb);

%% create some PTDF matrices
H   = makePTDF(baseMVA, bus, branch);
H1  = makePTDF(baseMVA, bus, branch, 1);
H4  = makePTDF(baseMVA, bus, branch, 4);
Heq = makePTDF(baseMVA, bus, branch, ones(nb, 1));
Hd  = makePTDF(baseMVA, bus, branch, Pd);
Hg  = makePTDF(baseMVA, bus, branch, Pg);

%% default slack
t_is(H, H1, 12, 'default slack');

%% using mpc
t_is(makePTDF(mpc), H, 12, 'MPC : default slack');
t_is(makePTDF(mpc, 1), H1, 12, 'H1 (from MPC)');
t_is(makePTDF(mpc, 4), H4, 12, 'H4 (from MPC)');
t_is(makePTDF(mpc, ones(nb,1)), Heq, 12, 'Heq (from MPC)');
t_is(makePTDF(mpc, Pd), Hd, 12, 'Hd (from MPC)');
t_is(makePTDF(mpc, Pg), Hg, 12, 'Hg (from MPC)');

%% vector slack same as scalar, using mpc
t_is(H1, makePTDF(mpc, e1), 12, 'H1 (scalar slack) = H1 (vector slack)');
t_is(H4, makePTDF(mpc, e4), 12, 'H4 (scalar slack) = H4 (vector slack)');

%% matrices get properly transformed by slack dist matrices
t_is(H1,  H1 * D1, 8,  'H1  == H1 * D1');
t_is(H4,  H1 * D4, 8,  'H4  == H1 * D4');
t_is(Heq, H1 * Deq, 8, 'Heq == H1 * Deq');
t_is(Hd,  H1 * Dd, 8,  'Hd  == H1 * Dd');
t_is(Hg,  H1 * Dg, 8,  'Hg  == H1 * Dg');
t_is(H1,  Heq * D1, 8,  'H1  == Heq * D1');
t_is(H4,  Heq * D4, 8,  'H4  == Heq * D4');
t_is(Heq, Heq * Deq, 8, 'Heq == Heq * Deq');
t_is(Hd,  Heq * Dd, 8,  'Hd  == Heq * Dd');
t_is(Hg,  Heq * Dg, 8,  'Hg  == Heq * Dg');
t_is(H1,  Hd * D1, 8,  'H1  == Hd * D1');
t_is(H4,  Hd * D4, 8,  'H4  == Hd * D4');
t_is(Heq, Hd * Deq, 8, 'Heq == Hd * Deq');
t_is(Hd,  Hd * Dd, 8,  'Hd  == Hd * Dd');
t_is(Hg,  Hd * Dg, 8,  'Hg  == Hd * Dg');

%% PTDFs can reconstruct flows
t_is(F,  H1 * P,  3,  'Flow == H1  * P');
t_is(F,  H4 * P,  3,  'Flow == H4  * P');
t_is(F,  Heq * P, 3,  'Flow == Heq * P');
t_is(F,  Hd * P,  3,  'Flow == Hd  * P');
t_is(F,  Hg * P,  3,  'Flow == Hg  * P');

%% other
t_is(F,  Hd * Pg,  3,  'Flow == Hd  * Pg');
t_is(F,  Hg * (-Pd),  3,  'Flow == Hg  * (-Pd)');
t_is(zeros(nbr,1),  Hd * (-Pd),  3,  'zeros == Hd  * (-Pd)');
t_is(zeros(nbr,1),  Hg * Pg,  3,  'zeros == Hg  * Pg');

%% single column, single slack
for k = 1:nb
    Hk = makePTDF(baseMVA, bus, branch, 1, k);
    t_is(Hk, H1(:, k), 12, sprintf('H1 : column %d', k));
end
for k = 1:nb
    Hk = makePTDF(mpc, 4, k);
    t_is(Hk, H4(:, k), 12, sprintf('H4 : column %d', k));
end

%% multiple columns, distributed slack
Hk = makePTDF(baseMVA, bus, branch, ones(nb, 1), (1:nb)');
t_is(Hk, Heq, 12, 'Heq : all columns');

Hk = makePTDF(mpc, Pd, find(Pd));
t_is(Hk, Hd(:, find(Pd)), 12, 'Hd : Pd columns');

Hk = makePTDF(baseMVA, bus, branch, Pg, find(Pg));
t_is(Hk, Hg(:, find(Pg)), 12, 'Hg : Pg columns');

%% specific transfers
for k = 1:nb
    txfr = zeros(nb, 1);    txfr(4) = -1;   txfr(k) = txfr(k) + 1;
    H = makePTDF(mpc, 4, txfr);
    t_is(H, H4(:, k), 12, sprintf('H4 (txfr) : column %d', k));
end
txfr = eye(9,9); txfr(1, :) = txfr(1, :) - 1;
H = makePTDF(mpc, 1, txfr);
t_is(H, H1, 12, sprintf('H1 (txfr) : full', k));
txfr = eye(9,9); txfr(4, :) = txfr(4, :) - 1;
H = makePTDF(mpc, 4, txfr);
t_is(H, H4, 12, sprintf('H4 (txfr) : full', k));

%% matrix of slacks (all cols)
if have_feature('matlab') && have_feature('matlab', 'vnum') < 9.001
    t_skip(2, 'MATLAB < 9.1 does not handle matrix ./ row-vector properly');
else
    Dm = [e4 e1 ones(nb, 1) Pd Pg e4 ones(nb, 1) Pd Pg];
    eHm = [H4(:,1) H1(:,2) Heq(:,3) Hd(:,4) Hg(:,5) H4(:,6) Heq(:,7) Hd(:,8) Hg(:,9)];
    Hm = makePTDF(mpc, Dm);
    t_is(Hm, eHm, 12, 'H (matrix slack) = H (vector slacks) - all cols');
    txfr = eye(nb, nb) - Dm ./ sum(Dm);
    Ht = makePTDF(mpc, [], txfr);
    t_is(Hm, Ht, 12, 'H (matrix slack) = H (equiv transfers) - all cols');
end
t_end;
