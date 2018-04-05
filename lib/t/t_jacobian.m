function t_jacobian(quiet)
%T_JACOBIAN  Numerical tests of partial derivative code.

%   MATPOWER
%   Copyright (c) 2004-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

t_begin(64, quiet);

casefile = 'case30';

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% run powerflow to get solved case
mpopt = mpoption('verbose', 0, 'out.all', 0);
mpc = loadcase(casefile);
[baseMVA, bus, gen, branch, success, et] = runpf(mpc, mpopt);

%% switch to internal bus numbering and build admittance matrices
[i2e, bus, gen, branch] = ext2int(bus, gen, branch);
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
Sbus = makeSbus(baseMVA, bus, gen);  %% net injected power in p.u.
Ybus_full   = full(Ybus);
Yf_full     = full(Yf);
Yt_full     = full(Yt);
Vm = bus(:, VM);
Va = bus(:, VA) * pi/180;
V = Vm .* exp(1j * Va);
Vr = real(V);
Vi = imag(V);
f = branch(:, F_BUS);       %% list of "from" buses
t = branch(:, T_BUS);       %% list of "to" buses
nl = length(f);
nb = length(V);
VV = V * ones(1, nb);
SS = Sbus * ones(1, nb);
pert = 1e-8;


%%-----  run tests for polar, then cartesian, coordinate voltages  -----
for vcart = 0:1
    %%-----  create perturbed voltages  -----
    if vcart        %% cartesian coordinate voltages (V1=Vr, V2=Vi)
        coord = 'cartesian';
        vv = {'r', 'i'};
        V1p = (Vr*ones(1,nb) + pert*eye(nb,nb)) + 1j * Vi * ones(1,nb);
        V2p = Vr*ones(1,nb) + 1j * (Vi*ones(1,nb) + pert*eye(nb,nb));
    else            %% polar coordinate voltages (V1=Va, V2=Vm)
        coord = 'polar';
        vv = {'a', 'm'};
        V1p = (Vm*ones(1,nb)) .* (exp(1j * (Va*ones(1,nb) + pert*eye(nb,nb))));
        V2p = (Vm*ones(1,nb) + pert*eye(nb,nb)) .* (exp(1j * Va) * ones(1,nb));
    end

    %%-----  check dImis_dV code  -----
    %% full matrices
    [dImis_dV1_full, dImis_dV2_full] = dImis_dV(Sbus, Ybus_full, V, vcart);

    %% sparse matrices
    [dImis_dV1, dImis_dV2] = dImis_dV(Sbus, Ybus, V, vcart);
    dImis_dV1_sp = full(dImis_dV1);
    dImis_dV2_sp = full(dImis_dV2);

    %% compute numerically to compare
    num_dImis_dV1 = full( ( Ybus * V1p - conj(SS./V1p) - Ybus*VV + conj(SS./VV)) / pert );
    num_dImis_dV2 = full( ( Ybus * V2p - conj(SS./V2p) - Ybus*VV + conj(SS./VV)) / pert );

    t_is(dImis_dV1_sp, num_dImis_dV1, 5, sprintf('%s - dImis_dV%s (sparse)', coord, vv{1}));
    t_is(dImis_dV2_sp, num_dImis_dV2, 5, sprintf('%s - dImis_dV%s (sparse)', coord, vv{2}));
    t_is(dImis_dV1_full, num_dImis_dV1, 5, sprintf('%s - dImis_dV%s (full)', coord, vv{1}));
    t_is(dImis_dV2_full, num_dImis_dV2, 5, sprintf('%s - dImis_dV%s (full)', coord, vv{2}));


    %%-----  check dSbus_dV code  -----
    %% full matrices
    [dSbus_dV1_full, dSbus_dV2_full] = dSbus_dV(Ybus_full, V, vcart);

    %% sparse matrices
    [dSbus_dV1, dSbus_dV2] = dSbus_dV(Ybus, V, vcart);
    dSbus_dV1_sp = full(dSbus_dV1);
    dSbus_dV2_sp = full(dSbus_dV2);

    %% compute numerically to compare
    num_dSbus_dV1 = full( (V1p .* conj(Ybus * V1p) - VV .* conj(Ybus * VV)) / pert );
    num_dSbus_dV2 = full( (V2p .* conj(Ybus * V2p) - VV .* conj(Ybus * VV)) / pert );

    t_is(dSbus_dV1_sp, num_dSbus_dV1, 5, sprintf('%s - dSbus_dV%s (sparse)', coord, vv{1}));
    t_is(dSbus_dV2_sp, num_dSbus_dV2, 5, sprintf('%s - dSbus_dV%s (sparse)', coord, vv{2}));
    t_is(dSbus_dV1_full, num_dSbus_dV1, 5, sprintf('%s - dSbus_dV%s (full)', coord, vv{1}));
    t_is(dSbus_dV2_full, num_dSbus_dV2, 5, sprintf('%s - dSbus_dV%s (full)', coord, vv{2}));

    %%-----  check dIbr_dV code  -----
    %% full matrices
    [dIf_dV1_full, dIf_dV2_full, dIt_dV1_full, dIt_dV2_full, If, It] = dIbr_dV(branch, Yf_full, Yt_full, V, vcart);

    %% sparse matrices
    [dIf_dV1, dIf_dV2, dIt_dV1, dIt_dV2, If, It] = dIbr_dV(branch, Yf, Yt, V, vcart);
    dIf_dV1_sp = full(dIf_dV1);
    dIf_dV2_sp = full(dIf_dV2);
    dIt_dV1_sp = full(dIt_dV1);
    dIt_dV2_sp = full(dIt_dV2);

    %% compute numerically to compare
    num_dIf_dV1 = full( (Yf * V1p - Yf * VV) / pert );
    num_dIf_dV2 = full( (Yf * V2p - Yf * VV) / pert );
    num_dIt_dV1 = full( (Yt * V1p - Yt * VV) / pert );
    num_dIt_dV2 = full( (Yt * V2p - Yt * VV) / pert );

    t_is(dIf_dV1_sp, num_dIf_dV1, 5, sprintf('%s - dIf_dV%s (sparse)', coord, vv{1}));
    t_is(dIf_dV2_sp, num_dIf_dV2, 5, sprintf('%s - dIf_dV%s (sparse)', coord, vv{2}));
    t_is(dIt_dV1_sp, num_dIt_dV1, 5, sprintf('%s - dIt_dV%s (sparse)', coord, vv{1}));
    t_is(dIt_dV2_sp, num_dIt_dV2, 5, sprintf('%s - dIt_dV%s (sparse)', coord, vv{2}));
    t_is(dIf_dV1_full, num_dIf_dV1, 5, sprintf('%s - dIf_dV%s (full)', coord, vv{1}));
    t_is(dIf_dV2_full, num_dIf_dV2, 5, sprintf('%s - dIf_dV%s (full)', coord, vv{2}));
    t_is(dIt_dV1_full, num_dIt_dV1, 5, sprintf('%s - dIt_dV%s (full)', coord, vv{1}));
    t_is(dIt_dV2_full, num_dIt_dV2, 5, sprintf('%s - dIt_dV%s (full)', coord, vv{2}));

    %%-----  check dSbr_dV code  -----
    %% full matrices
    [dSf_dV1_full, dSf_dV2_full, dSt_dV1_full, dSt_dV2_full, Sf, St] = dSbr_dV(branch, Yf_full, Yt_full, V, vcart);

    %% sparse matrices
    [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
    dSf_dV1_sp = full(dSf_dV1);
    dSf_dV2_sp = full(dSf_dV2);
    dSt_dV1_sp = full(dSt_dV1);
    dSt_dV2_sp = full(dSt_dV2);

    %% compute numerically to compare
    V1pf = V1p(f,:);
    V2pf = V2p(f,:);
    V1pt = V1p(t,:);
    V2pt = V2p(t,:);
    Sf2 = (V(f)*ones(1,nb)) .* conj(Yf * VV);
    St2 = (V(t)*ones(1,nb)) .* conj(Yt * VV);
    S1pf = V1pf .* conj(Yf * V1p);
    S2pf = V2pf .* conj(Yf * V2p);
    S1pt = V1pt .* conj(Yt * V1p);
    S2pt = V2pt .* conj(Yt * V2p);

    num_dSf_dV1 = full( (S1pf - Sf2) / pert );
    num_dSf_dV2 = full( (S2pf - Sf2) / pert );
    num_dSt_dV1 = full( (S1pt - St2) / pert );
    num_dSt_dV2 = full( (S2pt - St2) / pert );

    t_is(dSf_dV1_sp, num_dSf_dV1, 5, sprintf('%s - dSf_dV%s (sparse)', coord, vv{1}));
    t_is(dSf_dV2_sp, num_dSf_dV2, 5, sprintf('%s - dSf_dV%s (sparse)', coord, vv{2}));
    t_is(dSt_dV1_sp, num_dSt_dV1, 5, sprintf('%s - dSt_dV%s (sparse)', coord, vv{1}));
    t_is(dSt_dV2_sp, num_dSt_dV2, 5, sprintf('%s - dSt_dV%s (sparse)', coord, vv{2}));
    t_is(dSf_dV1_full, num_dSf_dV1, 5, sprintf('%s - dSf_dV%s (full)', coord, vv{1}));
    t_is(dSf_dV2_full, num_dSf_dV2, 5, sprintf('%s - dSf_dV%s (full)', coord, vv{2}));
    t_is(dSt_dV1_full, num_dSt_dV1, 5, sprintf('%s - dSt_dV%s (full)', coord, vv{1}));
    t_is(dSt_dV2_full, num_dSt_dV2, 5, sprintf('%s - dSt_dV%s (full)', coord, vv{2}));

    %%-----  check dAbr_dV code  -----
    %% full matrices
    [dAf_dV1_full, dAf_dV2_full, dAt_dV1_full, dAt_dV2_full] = ...
                            dAbr_dV(dSf_dV1_full, dSf_dV2_full, dSt_dV1_full, dSt_dV2_full, Sf, St);
    %% sparse matrices
    [dAf_dV1, dAf_dV2, dAt_dV1, dAt_dV2] = ...
                            dAbr_dV(dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St);
    dAf_dV1_sp = full(dAf_dV1);
    dAf_dV2_sp = full(dAf_dV2);
    dAt_dV1_sp = full(dAt_dV1);
    dAt_dV2_sp = full(dAt_dV2);

    %% compute numerically to compare
    num_dAf_dV1 = full( (abs(S1pf).^2 - abs(Sf2).^2) / pert );
    num_dAf_dV2 = full( (abs(S2pf).^2 - abs(Sf2).^2) / pert );
    num_dAt_dV1 = full( (abs(S1pt).^2 - abs(St2).^2) / pert );
    num_dAt_dV2 = full( (abs(S2pt).^2 - abs(St2).^2) / pert );

    t_is(dAf_dV1_sp, num_dAf_dV1, 4, sprintf('%s - dAf_dV%s (sparse)', coord, vv{1}));
    t_is(dAf_dV2_sp, num_dAf_dV2, 4, sprintf('%s - dAf_dV%s (sparse)', coord, vv{2}));
    t_is(dAt_dV1_sp, num_dAt_dV1, 4, sprintf('%s - dAt_dV%s (sparse)', coord, vv{1}));
    t_is(dAt_dV2_sp, num_dAt_dV2, 4, sprintf('%s - dAt_dV%s (sparse)', coord, vv{2}));
    t_is(dAf_dV1_full, num_dAf_dV1, 4, sprintf('%s - dAf_dV%s (full)', coord, vv{1}));
    t_is(dAf_dV2_full, num_dAf_dV2, 4, sprintf('%s - dAf_dV%s (full)', coord, vv{2}));
    t_is(dAt_dV1_full, num_dAt_dV1, 4, sprintf('%s - dAt_dV%s (full)', coord, vv{1}));
    t_is(dAt_dV2_full, num_dAt_dV2, 4, sprintf('%s - dAt_dV%s (full)', coord, vv{2}));
end

t_end;
