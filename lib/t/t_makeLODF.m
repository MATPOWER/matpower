function t_makeLODF(quiet)
% t_makeLODF - Tests for makeLODF.

%   MATPOWER
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

ntests = 35;
t_begin(ntests, quiet);

casefile = 't_auction_case';
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
mpc = loadcase(casefile);
mpopt = mpoption('out.all', 0, 'verbose', verbose);
[baseMVA, bus, gen, gencost, branch, f, success, et] = ...
    rundcopf(mpc, mpopt);
[i2e, bus, gen, branch] = ext2int(bus, gen, branch);

%% compute injections and flows
F0 = branch(:, PF);

%% create some PTDF matrices
H = makePTDF(baseMVA, bus, branch, 1);

%% create some PTDF matrices
s = warning('query', 'MATLAB:divideByZero');
warning('off', 'MATLAB:divideByZero');
LODF = makeLODF(branch, H);
warning(s.state, 'MATLAB:divideByZero');

%% take out non-essential lines one-by-one and see what happens
mpc.bus = bus;
mpc.gen = gen;
branch0 = branch;
outages = [1:12 14:15 17:18 20 27:33 35:41];
for k = outages
    mpc.branch = branch0;
    mpc.branch(k, BR_STATUS) = 0;
    [baseMVA, bus, gen, branch, success, et] = rundcpf(mpc, mpopt);
    F = branch(:, PF);

    t_is(LODF(:, k), (F - F0) / F0(k), 6, sprintf('LODF(:, %d)', k));
end

%% check with masking bridge branches
mpc.branch = branch0;
[~, bridges, nonbridges] = find_bridges(mpc);
t_is(length(bridges), 1, 12, 'length(bridges)')
t_is(bridges{1}, [13;16;19;21;22;23;24;34], 12, 'bridges');
t_is(nonbridges{1}, [1:12 14:15 17:18 20 27:33 35:41]', 12, 'nonbridges');
LODFm = makeLODF(mpc, H, 1);
LODF(:, bridges{1}) = NaN;
LODF([25:26], nonbridges{1}) = NaN;
t_is(LODFm, LODF, 12, 'with mask_bridge');

t_end;
