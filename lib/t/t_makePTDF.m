function t_makePTDF(quiet)
%T_MAKEPTDF  Tests for MAKEPTDF.

%   MATPOWER
%   Copyright (c) 2006-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

ntests = 24;
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
[baseMVA, bus, gen, gencost, branch, f, success, et] = ...
    rundcopf(casefile, mpopt);
[i2e, bus, gen, branch] = ext2int(bus, gen, branch);
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
Dg  = eye(nb) - Pd/sum(Pd) * ones(1, nb);
Dd  = eye(nb) - Pg/sum(Pg) * ones(1, nb);

%% create some PTDF matrices
H1  = makePTDF(baseMVA, bus, branch, 1);
H4  = makePTDF(baseMVA, bus, branch, 4);
Heq = makePTDF(baseMVA, bus, branch, ones(nb, 1));
Hg  = makePTDF(baseMVA, bus, branch, Pd);
Hd  = makePTDF(baseMVA, bus, branch, Pg);

%% matrices get properly transformed by slack dist matrices
t_is(H1,  H1 * D1, 8,  'H1  == H1 * D1');
t_is(H4,  H1 * D4, 8,  'H4  == H1 * D4');
t_is(Heq, H1 * Deq, 8, 'Heq == H1 * Deq');
t_is(Hg,  H1 * Dg, 8,  'Hg  == H1 * Dg');
t_is(Hd,  H1 * Dd, 8,  'Hd  == H1 * Dd');
t_is(H1,  Heq * D1, 8,  'H1  == Heq * D1');
t_is(H4,  Heq * D4, 8,  'H4  == Heq * D4');
t_is(Heq, Heq * Deq, 8, 'Heq == Heq * Deq');
t_is(Hg,  Heq * Dg, 8,  'Hg  == Heq * Dg');
t_is(Hd,  Heq * Dd, 8,  'Hd  == Heq * Dd');
t_is(H1,  Hg * D1, 8,  'H1  == Hg * D1');
t_is(H4,  Hg * D4, 8,  'H4  == Hg * D4');
t_is(Heq, Hg * Deq, 8, 'Heq == Hg * Deq');
t_is(Hg,  Hg * Dg, 8,  'Hg  == Hg * Dg');
t_is(Hd,  Hg * Dd, 8,  'Hd  == Hg * Dd');

%% PTDFs can reconstruct flows
t_is(F,  H1 * P,  3,  'Flow == H1  * P');
t_is(F,  H4 * P,  3,  'Flow == H4  * P');
t_is(F,  Heq * P, 3,  'Flow == Heq * P');
t_is(F,  Hg * P,  3,  'Flow == Hg  * P');
t_is(F,  Hd * P,  3,  'Flow == Hd  * P');

%% other
t_is(F,  Hg * Pg,  3,  'Flow == Hg  * Pg');
t_is(F,  Hd * (-Pd),  3,  'Flow == Hd  * (-Pd)');
t_is(zeros(nbr,1),  Hg * (-Pd),  3,  'zeros == Hg  * (-Pd)');
t_is(zeros(nbr,1),  Hd * Pg,  3,  'zeros == Hd  * Pg');

t_end;
