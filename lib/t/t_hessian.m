function t_hessian(quiet)
%T_HESSIAN  Numerical tests of 2nd derivative code.

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

t_begin(96, quiet);

casefile = 'case30';

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% run powerflow to get solved case
mpopt = mpoption('verbose', 0, 'out.all', 0);
[baseMVA, bus, gen, branch, success, et] = runpf(casefile, mpopt);

%% switch to internal bus numbering and build admittance matrices
[i2e, bus, gen, branch] = ext2int(bus, gen, branch);
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
Sbus = makeSbus(baseMVA, bus, gen);  %% net injected power in p.u.
Vm = bus(:, VM);
Va = bus(:, VA) * pi/180;
V = Vm .* exp(1j * Va);
Vr = real(V);
Vi = imag(V);
f = branch(:, F_BUS);       %% list of "from" buses
t = branch(:, T_BUS);       %% list of "to" buses
nl = length(f);
nb = length(V);
Cf = sparse(1:nl, f, ones(nl, 1), nl, nb);      %% connection matrix for line & from buses
Ct = sparse(1:nl, t, ones(nl, 1), nl, nb);      %% connection matrix for line & to buses
pert = 1e-8;

%%-----  run tests for polar, then cartesian, coordinate voltages  -----
for vcart = 0:1
    %%-----  create perturbed voltages  -----
    if vcart        %% cartesian coordinate voltages (V1=Vr, V2=Vi)
        coord = 'cartesian';
        vv = {'rr', 'ri', 'ir', 'ii'};
        V1p = (Vr*ones(1,nb) + pert*eye(nb,nb)) + 1j * Vi * ones(1,nb);
        V2p = Vr*ones(1,nb) + 1j * (Vi*ones(1,nb) + pert*eye(nb,nb));
    else            %% polar coordinate voltages (V1=Va, V2=Vm)
        coord = 'polar';
        vv = {'aa', 'av', 'va', 'vv'};
        V1p = (Vm*ones(1,nb)) .* (exp(1j * (Va*ones(1,nb) + pert*eye(nb,nb))));
        V2p = (Vm*ones(1,nb) + pert*eye(nb,nb)) .* (exp(1j * Va) * ones(1,nb));
    end


    %%-----  check d2Imis_dV2 code  -----
    t = ' - d2Imis_dV2 (complex current injections)';
    lam = 10 * rand(nb, 1);
    num_H11 = zeros(nb, nb);
    num_H12 = zeros(nb, nb);
    num_H21 = zeros(nb, nb);
    num_H22 = zeros(nb, nb);
    [dImis_dV1, dImis_dV2] = dImis_dV(Sbus, Ybus, V, vcart);
    [H11, H12, H21, H22] = d2Imis_dV2(Sbus, Ybus, V, lam, vcart);
    for i = 1:nb
        V1p = V;
        V2p = V;
        if vcart
            V1p(i) = (Vr(i) + pert) + 1j * Vi(i);       %% perturb Vr
            V2p(i) = Vr(i) + 1j * (Vi(i) + pert);       %% perturb Vi
        else
            V1p(i) = Vm(i) * exp(1j * (Va(i) + pert));  %% perturb Va
            V2p(i) = (Vm(i) + pert) * exp(1j * Va(i));  %% perturb Vm
        end
        [dImis_dV1_1p, dImis_dV2_1p] = dImis_dV(Sbus, Ybus, V1p, vcart);
        num_H11(:, i) = (dImis_dV1_1p - dImis_dV1).' * lam / pert;
        num_H21(:, i) = (dImis_dV2_1p - dImis_dV2).' * lam / pert;

        [dImis_dV1_2p, dImis_dV2_2p] = dImis_dV(Sbus, Ybus, V2p, vcart);
        num_H12(:, i) = (dImis_dV1_2p - dImis_dV1).' * lam / pert;
        num_H22(:, i) = (dImis_dV2_2p - dImis_dV2).' * lam / pert;
    end

    t_is(full(H11), num_H11, 4, sprintf('%s - H%s%s', coord, vv{1}, t));
    t_is(full(H12), num_H12, 4, sprintf('%s - H%s%s', coord, vv{2}, t));
    t_is(full(H21), num_H21, 4, sprintf('%s - H%s%s', coord, vv{3}, t));
    t_is(full(H22), num_H22, 4, sprintf('%s - H%s%s', coord, vv{4}, t));

    %%-----  check d2Sbus_dV2 code  -----
    t = ' - d2Sbus_dV2 (complex power injections)';
    lam = 10 * rand(nb, 1);
    num_H11 = zeros(nb, nb);
    num_H12 = zeros(nb, nb);
    num_H21 = zeros(nb, nb);
    num_H22 = zeros(nb, nb);
    [dSbus_dV1, dSbus_dV2] = dSbus_dV(Ybus, V, vcart);
    [H11, H12, H21, H22] = d2Sbus_dV2(Ybus, V, lam, vcart);
    for i = 1:nb
        V1p = V;
        V2p = V;
        if vcart
            V1p(i) = (Vr(i) + pert) + 1j * Vi(i);       %% perturb Vr
            V2p(i) = Vr(i) + 1j * (Vi(i) + pert);       %% perturb Vi
        else
            V1p(i) = Vm(i) * exp(1j * (Va(i) + pert));  %% perturb Va
            V2p(i) = (Vm(i) + pert) * exp(1j * Va(i));  %% perturb Vm
        end
        [dSbus_dV1_1p, dSbus_dV2_1p] = dSbus_dV(Ybus, V1p, vcart);
        num_H11(:, i) = (dSbus_dV1_1p - dSbus_dV1).' * lam / pert;
        num_H21(:, i) = (dSbus_dV2_1p - dSbus_dV2).' * lam / pert;

        [dSbus_dV1_2p, dSbus_dV2_2p] = dSbus_dV(Ybus, V2p, vcart);
        num_H12(:, i) = (dSbus_dV1_2p - dSbus_dV1).' * lam / pert;
        num_H22(:, i) = (dSbus_dV2_2p - dSbus_dV2).' * lam / pert;
    end

    t_is(full(H11), num_H11, 4, sprintf('%s - H%s%s', coord, vv{1}, t));
    t_is(full(H12), num_H12, 4, sprintf('%s - H%s%s', coord, vv{2}, t));
    t_is(full(H21), num_H21, 4, sprintf('%s - H%s%s', coord, vv{3}, t));
    t_is(full(H22), num_H22, 4, sprintf('%s - H%s%s', coord, vv{4}, t));

    %%-----  check d2Ibr_dV2 code  -----
    t = ' - d2Ibr_dV2 (complex currents)';
    lam = 10 * rand(nl, 1);
    % lam = [1; zeros(nl-1, 1)];
    num_Gf11 = zeros(nb, nb);
    num_Gf12 = zeros(nb, nb);
    num_Gf21 = zeros(nb, nb);
    num_Gf22 = zeros(nb, nb);
    num_Gt11 = zeros(nb, nb);
    num_Gt12 = zeros(nb, nb);
    num_Gt21 = zeros(nb, nb);
    num_Gt22 = zeros(nb, nb);
    [dIf_dV1, dIf_dV2, dIt_dV1, dIt_dV2, If, It] = dIbr_dV(branch, Yf, Yt, V, vcart);
    [Gf11, Gf12, Gf21, Gf22] = d2Ibr_dV2(Yf, V, lam, vcart);
    [Gt11, Gt12, Gt21, Gt22] = d2Ibr_dV2(Yt, V, lam, vcart);
    for i = 1:nb
        V1p = V;
        V2p = V;
        if vcart
            V1p(i) = (Vr(i) + pert) + 1j * Vi(i);       %% perturb Vr
            V2p(i) = Vr(i) + 1j * (Vi(i) + pert);       %% perturb Vi
        else
            V1p(i) = Vm(i) * exp(1j * (Va(i) + pert));  %% perturb Va
            V2p(i) = (Vm(i) + pert) * exp(1j * Va(i));  %% perturb Vm
        end
        [dIf_dV1_1p, dIf_dV2_1p, dIt_dV1_1p, dIt_dV2_1p, If_1p, It_1p] = ...
            dIbr_dV(branch, Yf, Yt, V1p, vcart);
        num_Gf11(:, i) = (dIf_dV1_1p - dIf_dV1).' * lam / pert;
        num_Gf21(:, i) = (dIf_dV2_1p - dIf_dV2).' * lam / pert;
        num_Gt11(:, i) = (dIt_dV1_1p - dIt_dV1).' * lam / pert;
        num_Gt21(:, i) = (dIt_dV2_1p - dIt_dV2).' * lam / pert;

        [dIf_dV1_2p, dIf_dV2_2p, dIt_dV1_2p, dIt_dV2_2p, If_2p, It_2p] = ...
            dIbr_dV(branch, Yf, Yt, V2p, vcart);
        num_Gf12(:, i) = (dIf_dV1_2p - dIf_dV1).' * lam / pert;
        num_Gf22(:, i) = (dIf_dV2_2p - dIf_dV2).' * lam / pert;
        num_Gt12(:, i) = (dIt_dV1_2p - dIt_dV1).' * lam / pert;
        num_Gt22(:, i) = (dIt_dV2_2p - dIt_dV2).' * lam / pert;
    end

    t_is(full(Gf11), num_Gf11, 4, sprintf('%s - Gf%s%s', coord, vv{1}, t));
    t_is(full(Gf12), num_Gf12, 4, sprintf('%s - Gf%s%s', coord, vv{2}, t));
    t_is(full(Gf21), num_Gf21, 4, sprintf('%s - Gf%s%s', coord, vv{3}, t));
    t_is(full(Gf22), num_Gf22, 4, sprintf('%s - Gf%s%s', coord, vv{4}, t));

    t_is(full(Gt11), num_Gt11, 4, sprintf('%s - Gt%s%s', coord, vv{1}, t));
    t_is(full(Gt12), num_Gt12, 4, sprintf('%s - Gt%s%s', coord, vv{2}, t));
    t_is(full(Gt21), num_Gt21, 4, sprintf('%s - Gt%s%s', coord, vv{3}, t));
    t_is(full(Gt22), num_Gt22, 4, sprintf('%s - Gt%s%s', coord, vv{4}, t));

    %%-----  check d2Sbr_dV2 code  -----
    t = ' - d2Sbr_dV2 (complex power flows)';
    lam = 10 * rand(nl, 1);
    % lam = [1; zeros(nl-1, 1)];
    num_Gf11 = zeros(nb, nb);
    num_Gf12 = zeros(nb, nb);
    num_Gf21 = zeros(nb, nb);
    num_Gf22 = zeros(nb, nb);
    num_Gt11 = zeros(nb, nb);
    num_Gt12 = zeros(nb, nb);
    num_Gt21 = zeros(nb, nb);
    num_Gt22 = zeros(nb, nb);
    [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
    [Gf11, Gf12, Gf21, Gf22] = d2Sbr_dV2(Cf, Yf, V, lam, vcart);
    [Gt11, Gt12, Gt21, Gt22] = d2Sbr_dV2(Ct, Yt, V, lam, vcart);
    for i = 1:nb
        V1p = V;
        V2p = V;
        if vcart
            V1p(i) = (Vr(i) + pert) + 1j * Vi(i);       %% perturb Vr
            V2p(i) = Vr(i) + 1j * (Vi(i) + pert);       %% perturb Vi
        else
            V1p(i) = Vm(i) * exp(1j * (Va(i) + pert));  %% perturb Va
            V2p(i) = (Vm(i) + pert) * exp(1j * Va(i));  %% perturb Vm
        end
        [dSf_dV1_1p, dSf_dV2_1p, dSt_dV1_1p, dSt_dV2_1p, Sf_1p, St_1p] = ...
            dSbr_dV(branch, Yf, Yt, V1p, vcart);
        num_Gf11(:, i) = (dSf_dV1_1p - dSf_dV1).' * lam / pert;
        num_Gf21(:, i) = (dSf_dV2_1p - dSf_dV2).' * lam / pert;
        num_Gt11(:, i) = (dSt_dV1_1p - dSt_dV1).' * lam / pert;
        num_Gt21(:, i) = (dSt_dV2_1p - dSt_dV2).' * lam / pert;

        [dSf_dV1_2p, dSf_dV2_2p, dSt_dV1_2p, dSt_dV2_2p, Sf_2p, St_2p] = ...
            dSbr_dV(branch, Yf, Yt, V2p, vcart);
        num_Gf12(:, i) = (dSf_dV1_2p - dSf_dV1).' * lam / pert;
        num_Gf22(:, i) = (dSf_dV2_2p - dSf_dV2).' * lam / pert;
        num_Gt12(:, i) = (dSt_dV1_2p - dSt_dV1).' * lam / pert;
        num_Gt22(:, i) = (dSt_dV2_2p - dSt_dV2).' * lam / pert;
    end

    t_is(full(Gf11), num_Gf11, 4, sprintf('%s - Gf%s%s', coord, vv{1}, t));
    t_is(full(Gf12), num_Gf12, 4, sprintf('%s - Gf%s%s', coord, vv{2}, t));
    t_is(full(Gf21), num_Gf21, 4, sprintf('%s - Gf%s%s', coord, vv{3}, t));
    t_is(full(Gf22), num_Gf22, 4, sprintf('%s - Gf%s%s', coord, vv{4}, t));

    t_is(full(Gt11), num_Gt11, 4, sprintf('%s - Gt%s%s', coord, vv{1}, t));
    t_is(full(Gt12), num_Gt12, 4, sprintf('%s - Gt%s%s', coord, vv{2}, t));
    t_is(full(Gt21), num_Gt21, 4, sprintf('%s - Gt%s%s', coord, vv{3}, t));
    t_is(full(Gt22), num_Gt22, 4, sprintf('%s - Gt%s%s', coord, vv{4}, t));

    %%-----  check d2Abr_dV2 code  -----
    t = ' - d2Abr_dV2 (squared current magnitudes)';
    lam = 10 * rand(nl, 1);
    % lam = [1; zeros(nl-1, 1)];
    num_Gf11 = zeros(nb, nb);
    num_Gf12 = zeros(nb, nb);
    num_Gf21 = zeros(nb, nb);
    num_Gf22 = zeros(nb, nb);
    num_Gt11 = zeros(nb, nb);
    num_Gt12 = zeros(nb, nb);
    num_Gt21 = zeros(nb, nb);
    num_Gt22 = zeros(nb, nb);
    d2If_dV2 = @(V, mu)d2Ibr_dV2(Yf, V, mu, vcart);
    d2It_dV2 = @(V, mu)d2Ibr_dV2(Yt, V, mu, vcart);
    [dIf_dV1, dIf_dV2, dIt_dV1, dIt_dV2, If, It] = dIbr_dV(branch, Yf, Yt, V, vcart);
    [dAf_dV1, dAf_dV2, dAt_dV1, dAt_dV2] = ...
                            dAbr_dV(dIf_dV1, dIf_dV2, dIt_dV1, dIt_dV2, If, It);
    [Gf11, Gf12, Gf21, Gf22] = d2Abr_dV2(d2If_dV2, dIf_dV1, dIf_dV2, If, V, lam);
    [Gt11, Gt12, Gt21, Gt22] = d2Abr_dV2(d2It_dV2, dIt_dV1, dIt_dV2, It, V, lam);
    for i = 1:nb
        V1p = V;
        V2p = V;
        if vcart
            V1p(i) = (Vr(i) + pert) + 1j * Vi(i);       %% perturb Vr
            V2p(i) = Vr(i) + 1j * (Vi(i) + pert);       %% perturb Vi
        else
            V1p(i) = Vm(i) * exp(1j * (Va(i) + pert));  %% perturb Va
            V2p(i) = (Vm(i) + pert) * exp(1j * Va(i));  %% perturb Vm
        end
        [dIf_dV1_1p, dIf_dV2_1p, dIt_dV1_1p, dIt_dV2_1p, If_1p, It_1p] = ...
            dIbr_dV(branch, Yf, Yt, V1p, vcart);
        [dAf_dV1_1p, dAf_dV2_1p, dAt_dV1_1p, dAt_dV2_1p] = ...
            dAbr_dV(dIf_dV1_1p, dIf_dV2_1p, dIt_dV1_1p, dIt_dV2_1p, If_1p, It_1p);
        num_Gf11(:, i) = (dAf_dV1_1p - dAf_dV1).' * lam / pert;
        num_Gf21(:, i) = (dAf_dV2_1p - dAf_dV2).' * lam / pert;
        num_Gt11(:, i) = (dAt_dV1_1p - dAt_dV1).' * lam / pert;
        num_Gt21(:, i) = (dAt_dV2_1p - dAt_dV2).' * lam / pert;

        [dIf_dV1_2p, dIf_dV2_2p, dIt_dV1_2p, dIt_dV2_2p, If_2p, It_2p] = ...
            dIbr_dV(branch, Yf, Yt, V2p, vcart);
        [dAf_dV1_2p, dAf_dV2_2p, dAt_dV1_2p, dAt_dV2_2p] = ...
            dAbr_dV(dIf_dV1_2p, dIf_dV2_2p, dIt_dV1_2p, dIt_dV2_2p, If_2p, It_2p);
        num_Gf12(:, i) = (dAf_dV1_2p - dAf_dV1).' * lam / pert;
        num_Gf22(:, i) = (dAf_dV2_2p - dAf_dV2).' * lam / pert;
        num_Gt12(:, i) = (dAt_dV1_2p - dAt_dV1).' * lam / pert;
        num_Gt22(:, i) = (dAt_dV2_2p - dAt_dV2).' * lam / pert;
    end

    t_is(full(Gf11), num_Gf11, 3, sprintf('%s - Gf%s%s', coord, vv{1}, t));
    t_is(full(Gf12), num_Gf12, 3, sprintf('%s - Gf%s%s', coord, vv{2}, t));
    t_is(full(Gf21), num_Gf21, 3, sprintf('%s - Gf%s%s', coord, vv{3}, t));
    t_is(full(Gf22), num_Gf22, 3, sprintf('%s - Gf%s%s', coord, vv{4}, t));

    t_is(full(Gt11), num_Gt11, 3, sprintf('%s - Gt%s%s', coord, vv{1}, t));
    t_is(full(Gt12), num_Gt12, 3, sprintf('%s - Gt%s%s', coord, vv{2}, t));
    t_is(full(Gt21), num_Gt21, 3, sprintf('%s - Gt%s%s', coord, vv{3}, t));
    t_is(full(Gt22), num_Gt22, 3, sprintf('%s - Gt%s%s', coord, vv{4}, t));

    %%-----  check d2Abr_dV2 code  -----
    t = ' - d2Abr_dV2 (squared apparent power flows)';
    lam = 10 * rand(nl, 1);
    % lam = [1; zeros(nl-1, 1)];
    num_Gf11 = zeros(nb, nb);
    num_Gf12 = zeros(nb, nb);
    num_Gf21 = zeros(nb, nb);
    num_Gf22 = zeros(nb, nb);
    num_Gt11 = zeros(nb, nb);
    num_Gt12 = zeros(nb, nb);
    num_Gt21 = zeros(nb, nb);
    num_Gt22 = zeros(nb, nb);
    d2Sf_dV2 = @(V, mu)d2Sbr_dV2(Cf, Yf, V, mu, vcart);
    d2St_dV2 = @(V, mu)d2Sbr_dV2(Ct, Yt, V, mu, vcart);
    [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
    [dAf_dV1, dAf_dV2, dAt_dV1, dAt_dV2] = ...
                            dAbr_dV(dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St);
    [Gf11, Gf12, Gf21, Gf22] = d2Abr_dV2(d2Sf_dV2, dSf_dV1, dSf_dV2, Sf, V, lam);
    [Gt11, Gt12, Gt21, Gt22] = d2Abr_dV2(d2St_dV2, dSt_dV1, dSt_dV2, St, V, lam);
    for i = 1:nb
        V1p = V;
        V2p = V;
        if vcart
            V1p(i) = (Vr(i) + pert) + 1j * Vi(i);       %% perturb Vr
            V2p(i) = Vr(i) + 1j * (Vi(i) + pert);       %% perturb Vi
        else
            V1p(i) = Vm(i) * exp(1j * (Va(i) + pert));  %% perturb Va
            V2p(i) = (Vm(i) + pert) * exp(1j * Va(i));  %% perturb Vm
        end
        [dSf_dV1_1p, dSf_dV2_1p, dSt_dV1_1p, dSt_dV2_1p, Sf_1p, St_1p] = ...
            dSbr_dV(branch, Yf, Yt, V1p, vcart);
        [dAf_dV1_1p, dAf_dV2_1p, dAt_dV1_1p, dAt_dV2_1p] = ...
            dAbr_dV(dSf_dV1_1p, dSf_dV2_1p, dSt_dV1_1p, dSt_dV2_1p, Sf_1p, St_1p);
        num_Gf11(:, i) = (dAf_dV1_1p - dAf_dV1).' * lam / pert;
        num_Gf21(:, i) = (dAf_dV2_1p - dAf_dV2).' * lam / pert;
        num_Gt11(:, i) = (dAt_dV1_1p - dAt_dV1).' * lam / pert;
        num_Gt21(:, i) = (dAt_dV2_1p - dAt_dV2).' * lam / pert;

        [dSf_dV1_2p, dSf_dV2_2p, dSt_dV1_2p, dSt_dV2_2p, Sf_2p, St_2p] = ...
            dSbr_dV(branch, Yf, Yt, V2p, vcart);
        [dAf_dV1_2p, dAf_dV2_2p, dAt_dV1_2p, dAt_dV2_2p] = ...
            dAbr_dV(dSf_dV1_2p, dSf_dV2_2p, dSt_dV1_2p, dSt_dV2_2p, Sf_2p, St_2p);
        num_Gf12(:, i) = (dAf_dV1_2p - dAf_dV1).' * lam / pert;
        num_Gf22(:, i) = (dAf_dV2_2p - dAf_dV2).' * lam / pert;
        num_Gt12(:, i) = (dAt_dV1_2p - dAt_dV1).' * lam / pert;
        num_Gt22(:, i) = (dAt_dV2_2p - dAt_dV2).' * lam / pert;
    end

    t_is(full(Gf11), num_Gf11, 2, sprintf('%s - Gf%s%s', coord, vv{1}, t));
    t_is(full(Gf12), num_Gf12, 2, sprintf('%s - Gf%s%s', coord, vv{2}, t));
    t_is(full(Gf21), num_Gf21, 2, sprintf('%s - Gf%s%s', coord, vv{3}, t));
    t_is(full(Gf22), num_Gf22, 2, sprintf('%s - Gf%s%s', coord, vv{4}, t));

    t_is(full(Gt11), num_Gt11, 2, sprintf('%s - Gt%s%s', coord, vv{1}, t));
    t_is(full(Gt12), num_Gt12, 2, sprintf('%s - Gt%s%s', coord, vv{2}, t));
    t_is(full(Gt21), num_Gt21, 2, sprintf('%s - Gt%s%s', coord, vv{3}, t));
    t_is(full(Gt22), num_Gt22, 2, sprintf('%s - Gt%s%s', coord, vv{4}, t));

    %%-----  check d2Abr_dV2 code  -----
    t = ' - d2Abr_dV2 (squared real power flows)';
    lam = 10 * rand(nl, 1);
    % lam = [1; zeros(nl-1, 1)];
    num_Gf11 = zeros(nb, nb);
    num_Gf12 = zeros(nb, nb);
    num_Gf21 = zeros(nb, nb);
    num_Gf22 = zeros(nb, nb);
    num_Gt11 = zeros(nb, nb);
    num_Gt12 = zeros(nb, nb);
    num_Gt21 = zeros(nb, nb);
    num_Gt22 = zeros(nb, nb);
    d2Sf_dV2 = @(V, mu)d2Sbr_dV2(Cf, Yf, V, mu, vcart);
    d2St_dV2 = @(V, mu)d2Sbr_dV2(Ct, Yt, V, mu, vcart);
    [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
    [dAf_dV1, dAf_dV2, dAt_dV1, dAt_dV2] = ...
                            dAbr_dV(real(dSf_dV1), real(dSf_dV2), real(dSt_dV1), real(dSt_dV2), real(Sf), real(St));
    [Gf11, Gf12, Gf21, Gf22] = d2Abr_dV2(d2Sf_dV2, real(dSf_dV1), real(dSf_dV2), real(Sf), V, lam);
    [Gt11, Gt12, Gt21, Gt22] = d2Abr_dV2(d2St_dV2, real(dSt_dV1), real(dSt_dV2), real(St), V, lam);
    for i = 1:nb
        V1p = V;
        V2p = V;
        if vcart
            V1p(i) = (Vr(i) + pert) + 1j * Vi(i);       %% perturb Vr
            V2p(i) = Vr(i) + 1j * (Vi(i) + pert);       %% perturb Vi
        else
            V1p(i) = Vm(i) * exp(1j * (Va(i) + pert));  %% perturb Va
            V2p(i) = (Vm(i) + pert) * exp(1j * Va(i));  %% perturb Vm
        end
        [dSf_dV1_1p, dSf_dV2_1p, dSt_dV1_1p, dSt_dV2_1p, Sf_1p, St_1p] = ...
            dSbr_dV(branch, Yf, Yt, V1p, vcart);
        [dAf_dV1_1p, dAf_dV2_1p, dAt_dV1_1p, dAt_dV2_1p] = ...
            dAbr_dV(real(dSf_dV1_1p), real(dSf_dV2_1p), real(dSt_dV1_1p), real(dSt_dV2_1p), real(Sf_1p), real(St_1p));
        num_Gf11(:, i) = (dAf_dV1_1p - dAf_dV1).' * lam / pert;
        num_Gf21(:, i) = (dAf_dV2_1p - dAf_dV2).' * lam / pert;
        num_Gt11(:, i) = (dAt_dV1_1p - dAt_dV1).' * lam / pert;
        num_Gt21(:, i) = (dAt_dV2_1p - dAt_dV2).' * lam / pert;

        [dSf_dV1_2p, dSf_dV2_2p, dSt_dV1_2p, dSt_dV2_2p, Sf_2p, St_2p] = ...
            dSbr_dV(branch, Yf, Yt, V2p, vcart);
        [dAf_dV1_2p, dAf_dV2_2p, dAt_dV1_2p, dAt_dV2_2p] = ...
            dAbr_dV(real(dSf_dV1_2p), real(dSf_dV2_2p), real(dSt_dV1_2p), real(dSt_dV2_2p), real(Sf_2p), real(St_2p));
        num_Gf12(:, i) = (dAf_dV1_2p - dAf_dV1).' * lam / pert;
        num_Gf22(:, i) = (dAf_dV2_2p - dAf_dV2).' * lam / pert;
        num_Gt12(:, i) = (dAt_dV1_2p - dAt_dV1).' * lam / pert;
        num_Gt22(:, i) = (dAt_dV2_2p - dAt_dV2).' * lam / pert;
    end

    t_is(full(Gf11), num_Gf11, 2, sprintf('%s - Gf%s%s', coord, vv{1}, t));
    t_is(full(Gf12), num_Gf12, 2, sprintf('%s - Gf%s%s', coord, vv{2}, t));
    t_is(full(Gf21), num_Gf21, 2, sprintf('%s - Gf%s%s', coord, vv{3}, t));
    t_is(full(Gf22), num_Gf22, 2, sprintf('%s - Gf%s%s', coord, vv{4}, t));

    t_is(full(Gt11), num_Gt11, 2, sprintf('%s - Gt%s%s', coord, vv{1}, t));
    t_is(full(Gt12), num_Gt12, 2, sprintf('%s - Gt%s%s', coord, vv{2}, t));
    t_is(full(Gt21), num_Gt21, 2, sprintf('%s - Gt%s%s', coord, vv{3}, t));
    t_is(full(Gt22), num_Gt22, 2, sprintf('%s - Gt%s%s', coord, vv{4}, t));
end

t_end;
