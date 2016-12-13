function t_hessian(quiet)
%T_HESSIAN  Numerical tests of 2nd derivative code.

%   MATPOWER
%   Copyright (c) 2004-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

t_begin(44, quiet);

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
Vm = bus(:, VM);
Va = bus(:, VA) * pi/180;
V = Vm .* exp(1j * Va);
f = branch(:, F_BUS);       %% list of "from" buses
t = branch(:, T_BUS);       %% list of "to" buses
nl = length(f);
nb = length(V);
Cf = sparse(1:nl, f, ones(nl, 1), nl, nb);      %% connection matrix for line & from buses
Ct = sparse(1:nl, t, ones(nl, 1), nl, nb);      %% connection matrix for line & to buses
pert = 1e-8;

%%-----  check d2Sbus_dV2 code  -----
t = ' - d2Sbus_dV2 (complex power injections)';
lam = 10 * rand(nb, 1);
num_Haa = zeros(nb, nb);
num_Hav = zeros(nb, nb);
num_Hva = zeros(nb, nb);
num_Hvv = zeros(nb, nb);
[dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V);
[Haa, Hav, Hva, Hvv] = d2Sbus_dV2(Ybus, V, lam);
for i = 1:nb
    Vap = V;
    Vap(i) = Vm(i) * exp(1j * (Va(i) + pert));
    [dSbus_dVm_ap, dSbus_dVa_ap] = dSbus_dV(Ybus, Vap);
    num_Haa(:, i) = (dSbus_dVa_ap - dSbus_dVa).' * lam / pert;
    num_Hva(:, i) = (dSbus_dVm_ap - dSbus_dVm).' * lam / pert;

    Vmp = V;
    Vmp(i) = (Vm(i) + pert) * exp(1j * Va(i));
    [dSbus_dVm_mp, dSbus_dVa_mp] = dSbus_dV(Ybus, Vmp);
    num_Hav(:, i) = (dSbus_dVa_mp - dSbus_dVa).' * lam / pert;
    num_Hvv(:, i) = (dSbus_dVm_mp - dSbus_dVm).' * lam / pert;
end

t_is(full(Haa), num_Haa, 4, ['Haa' t]);
t_is(full(Hav), num_Hav, 4, ['Hav' t]);
t_is(full(Hva), num_Hva, 4, ['Hva' t]);
t_is(full(Hvv), num_Hvv, 4, ['Hvv' t]);

%%-----  check d2Sbr_dV2 code  -----
t = ' - d2Sbr_dV2 (complex power flows)';
lam = 10 * rand(nl, 1);
% lam = [1; zeros(nl-1, 1)];
num_Gfaa = zeros(nb, nb);
num_Gfav = zeros(nb, nb);
num_Gfva = zeros(nb, nb);
num_Gfvv = zeros(nb, nb);
num_Gtaa = zeros(nb, nb);
num_Gtav = zeros(nb, nb);
num_Gtva = zeros(nb, nb);
num_Gtvv = zeros(nb, nb);
[dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = dSbr_dV(branch, Yf, Yt, V);
[Gfaa, Gfav, Gfva, Gfvv] = d2Sbr_dV2(Cf, Yf, V, lam);
[Gtaa, Gtav, Gtva, Gtvv] = d2Sbr_dV2(Ct, Yt, V, lam);
for i = 1:nb
    Vap = V;
    Vap(i) = Vm(i) * exp(1j * (Va(i) + pert));
    [dSf_dVa_ap, dSf_dVm_ap, dSt_dVa_ap, dSt_dVm_ap, Sf_ap, St_ap] = ...
        dSbr_dV(branch, Yf, Yt, Vap);
    num_Gfaa(:, i) = (dSf_dVa_ap - dSf_dVa).' * lam / pert;
    num_Gfva(:, i) = (dSf_dVm_ap - dSf_dVm).' * lam / pert;
    num_Gtaa(:, i) = (dSt_dVa_ap - dSt_dVa).' * lam / pert;
    num_Gtva(:, i) = (dSt_dVm_ap - dSt_dVm).' * lam / pert;

    Vmp = V;
    Vmp(i) = (Vm(i) + pert) * exp(1j * Va(i));
    [dSf_dVa_mp, dSf_dVm_mp, dSt_dVa_mp, dSt_dVm_mp, Sf_mp, St_mp] = ...
        dSbr_dV(branch, Yf, Yt, Vmp);
    num_Gfav(:, i) = (dSf_dVa_mp - dSf_dVa).' * lam / pert;
    num_Gfvv(:, i) = (dSf_dVm_mp - dSf_dVm).' * lam / pert;
    num_Gtav(:, i) = (dSt_dVa_mp - dSt_dVa).' * lam / pert;
    num_Gtvv(:, i) = (dSt_dVm_mp - dSt_dVm).' * lam / pert;
end

t_is(full(Gfaa), num_Gfaa, 4, ['Gfaa' t]);
t_is(full(Gfav), num_Gfav, 4, ['Gfav' t]);
t_is(full(Gfva), num_Gfva, 4, ['Gfva' t]);
t_is(full(Gfvv), num_Gfvv, 4, ['Gfvv' t]);

t_is(full(Gtaa), num_Gtaa, 4, ['Gtaa' t]);
t_is(full(Gtav), num_Gtav, 4, ['Gtav' t]);
t_is(full(Gtva), num_Gtva, 4, ['Gtva' t]);
t_is(full(Gtvv), num_Gtvv, 4, ['Gtvv' t]);

%%-----  check d2Ibr_dV2 code  -----
t = ' - d2Ibr_dV2 (complex currents)';
lam = 10 * rand(nl, 1);
% lam = [1; zeros(nl-1, 1)];
num_Gfaa = zeros(nb, nb);
num_Gfav = zeros(nb, nb);
num_Gfva = zeros(nb, nb);
num_Gfvv = zeros(nb, nb);
num_Gtaa = zeros(nb, nb);
num_Gtav = zeros(nb, nb);
num_Gtva = zeros(nb, nb);
num_Gtvv = zeros(nb, nb);
[dIf_dVa, dIf_dVm, dIt_dVa, dIt_dVm, If, It] = dIbr_dV(branch, Yf, Yt, V);
[Gfaa, Gfav, Gfva, Gfvv] = d2Ibr_dV2(Yf, V, lam);
[Gtaa, Gtav, Gtva, Gtvv] = d2Ibr_dV2(Yt, V, lam);
for i = 1:nb
    Vap = V;
    Vap(i) = Vm(i) * exp(1j * (Va(i) + pert));
    [dIf_dVa_ap, dIf_dVm_ap, dIt_dVa_ap, dIt_dVm_ap, If_ap, It_ap] = ...
        dIbr_dV(branch, Yf, Yt, Vap);
    num_Gfaa(:, i) = (dIf_dVa_ap - dIf_dVa).' * lam / pert;
    num_Gfva(:, i) = (dIf_dVm_ap - dIf_dVm).' * lam / pert;
    num_Gtaa(:, i) = (dIt_dVa_ap - dIt_dVa).' * lam / pert;
    num_Gtva(:, i) = (dIt_dVm_ap - dIt_dVm).' * lam / pert;

    Vmp = V;
    Vmp(i) = (Vm(i) + pert) * exp(1j * Va(i));
    [dIf_dVa_mp, dIf_dVm_mp, dIt_dVa_mp, dIt_dVm_mp, If_mp, It_mp] = ...
        dIbr_dV(branch, Yf, Yt, Vmp);
    num_Gfav(:, i) = (dIf_dVa_mp - dIf_dVa).' * lam / pert;
    num_Gfvv(:, i) = (dIf_dVm_mp - dIf_dVm).' * lam / pert;
    num_Gtav(:, i) = (dIt_dVa_mp - dIt_dVa).' * lam / pert;
    num_Gtvv(:, i) = (dIt_dVm_mp - dIt_dVm).' * lam / pert;
end

t_is(full(Gfaa), num_Gfaa, 4, ['Gfaa' t]);
t_is(full(Gfav), num_Gfav, 4, ['Gfav' t]);
t_is(full(Gfva), num_Gfva, 4, ['Gfva' t]);
t_is(full(Gfvv), num_Gfvv, 4, ['Gfvv' t]);

t_is(full(Gtaa), num_Gtaa, 4, ['Gtaa' t]);
t_is(full(Gtav), num_Gtav, 4, ['Gtav' t]);
t_is(full(Gtva), num_Gtva, 4, ['Gtva' t]);
t_is(full(Gtvv), num_Gtvv, 4, ['Gtvv' t]);

%%-----  check d2ASbr_dV2 code  -----
t = ' - d2ASbr_dV2 (squared apparent power flows)';
lam = 10 * rand(nl, 1);
% lam = [1; zeros(nl-1, 1)];
num_Gfaa = zeros(nb, nb);
num_Gfav = zeros(nb, nb);
num_Gfva = zeros(nb, nb);
num_Gfvv = zeros(nb, nb);
num_Gtaa = zeros(nb, nb);
num_Gtav = zeros(nb, nb);
num_Gtva = zeros(nb, nb);
num_Gtvv = zeros(nb, nb);
[dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = dSbr_dV(branch, Yf, Yt, V);
[dAf_dVa, dAf_dVm, dAt_dVa, dAt_dVm] = ...
                        dAbr_dV(dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St);
[Gfaa, Gfav, Gfva, Gfvv] = d2ASbr_dV2(dSf_dVa, dSf_dVm, Sf, Cf, Yf, V, lam);
[Gtaa, Gtav, Gtva, Gtvv] = d2ASbr_dV2(dSt_dVa, dSt_dVm, St, Ct, Yt, V, lam);
for i = 1:nb
    Vap = V;
    Vap(i) = Vm(i) * exp(1j * (Va(i) + pert));
    [dSf_dVa_ap, dSf_dVm_ap, dSt_dVa_ap, dSt_dVm_ap, Sf_ap, St_ap] = ...
        dSbr_dV(branch, Yf, Yt, Vap);
    [dAf_dVa_ap, dAf_dVm_ap, dAt_dVa_ap, dAt_dVm_ap] = ...
        dAbr_dV(dSf_dVa_ap, dSf_dVm_ap, dSt_dVa_ap, dSt_dVm_ap, Sf_ap, St_ap);
    num_Gfaa(:, i) = (dAf_dVa_ap - dAf_dVa).' * lam / pert;
    num_Gfva(:, i) = (dAf_dVm_ap - dAf_dVm).' * lam / pert;
    num_Gtaa(:, i) = (dAt_dVa_ap - dAt_dVa).' * lam / pert;
    num_Gtva(:, i) = (dAt_dVm_ap - dAt_dVm).' * lam / pert;

    Vmp = V;
    Vmp(i) = (Vm(i) + pert) * exp(1j * Va(i));
    [dSf_dVa_mp, dSf_dVm_mp, dSt_dVa_mp, dSt_dVm_mp, Sf_mp, St_mp] = ...
        dSbr_dV(branch, Yf, Yt, Vmp);
    [dAf_dVa_mp, dAf_dVm_mp, dAt_dVa_mp, dAt_dVm_mp] = ...
        dAbr_dV(dSf_dVa_mp, dSf_dVm_mp, dSt_dVa_mp, dSt_dVm_mp, Sf_mp, St_mp);
    num_Gfav(:, i) = (dAf_dVa_mp - dAf_dVa).' * lam / pert;
    num_Gfvv(:, i) = (dAf_dVm_mp - dAf_dVm).' * lam / pert;
    num_Gtav(:, i) = (dAt_dVa_mp - dAt_dVa).' * lam / pert;
    num_Gtvv(:, i) = (dAt_dVm_mp - dAt_dVm).' * lam / pert;
end

t_is(full(Gfaa), num_Gfaa, 2, ['Gfaa' t]);
t_is(full(Gfav), num_Gfav, 2, ['Gfav' t]);
t_is(full(Gfva), num_Gfva, 2, ['Gfva' t]);
t_is(full(Gfvv), num_Gfvv, 2, ['Gfvv' t]);

t_is(full(Gtaa), num_Gtaa, 2, ['Gtaa' t]);
t_is(full(Gtav), num_Gtav, 2, ['Gtav' t]);
t_is(full(Gtva), num_Gtva, 2, ['Gtva' t]);
t_is(full(Gtvv), num_Gtvv, 2, ['Gtvv' t]);

%%-----  check d2ASbr_dV2 code  -----
t = ' - d2ASbr_dV2 (squared real power flows)';
lam = 10 * rand(nl, 1);
% lam = [1; zeros(nl-1, 1)];
num_Gfaa = zeros(nb, nb);
num_Gfav = zeros(nb, nb);
num_Gfva = zeros(nb, nb);
num_Gfvv = zeros(nb, nb);
num_Gtaa = zeros(nb, nb);
num_Gtav = zeros(nb, nb);
num_Gtva = zeros(nb, nb);
num_Gtvv = zeros(nb, nb);
[dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = dSbr_dV(branch, Yf, Yt, V);
[dAf_dVa, dAf_dVm, dAt_dVa, dAt_dVm] = ...
                        dAbr_dV(real(dSf_dVa), real(dSf_dVm), real(dSt_dVa), real(dSt_dVm), real(Sf), real(St));
[Gfaa, Gfav, Gfva, Gfvv] = d2ASbr_dV2(real(dSf_dVa), real(dSf_dVm), real(Sf), Cf, Yf, V, lam);
[Gtaa, Gtav, Gtva, Gtvv] = d2ASbr_dV2(real(dSt_dVa), real(dSt_dVm), real(St), Ct, Yt, V, lam);
for i = 1:nb
    Vap = V;
    Vap(i) = Vm(i) * exp(1j * (Va(i) + pert));
    [dSf_dVa_ap, dSf_dVm_ap, dSt_dVa_ap, dSt_dVm_ap, Sf_ap, St_ap] = ...
        dSbr_dV(branch, Yf, Yt, Vap);
    [dAf_dVa_ap, dAf_dVm_ap, dAt_dVa_ap, dAt_dVm_ap] = ...
        dAbr_dV(real(dSf_dVa_ap), real(dSf_dVm_ap), real(dSt_dVa_ap), real(dSt_dVm_ap), real(Sf_ap), real(St_ap));
    num_Gfaa(:, i) = (dAf_dVa_ap - dAf_dVa).' * lam / pert;
    num_Gfva(:, i) = (dAf_dVm_ap - dAf_dVm).' * lam / pert;
    num_Gtaa(:, i) = (dAt_dVa_ap - dAt_dVa).' * lam / pert;
    num_Gtva(:, i) = (dAt_dVm_ap - dAt_dVm).' * lam / pert;

    Vmp = V;
    Vmp(i) = (Vm(i) + pert) * exp(1j * Va(i));
    [dSf_dVa_mp, dSf_dVm_mp, dSt_dVa_mp, dSt_dVm_mp, Sf_mp, St_mp] = ...
        dSbr_dV(branch, Yf, Yt, Vmp);
    [dAf_dVa_mp, dAf_dVm_mp, dAt_dVa_mp, dAt_dVm_mp] = ...
        dAbr_dV(real(dSf_dVa_mp), real(dSf_dVm_mp), real(dSt_dVa_mp), real(dSt_dVm_mp), real(Sf_mp), real(St_mp));
    num_Gfav(:, i) = (dAf_dVa_mp - dAf_dVa).' * lam / pert;
    num_Gfvv(:, i) = (dAf_dVm_mp - dAf_dVm).' * lam / pert;
    num_Gtav(:, i) = (dAt_dVa_mp - dAt_dVa).' * lam / pert;
    num_Gtvv(:, i) = (dAt_dVm_mp - dAt_dVm).' * lam / pert;
end

t_is(full(Gfaa), num_Gfaa, 2, ['Gfaa' t]);
t_is(full(Gfav), num_Gfav, 2, ['Gfav' t]);
t_is(full(Gfva), num_Gfva, 2, ['Gfva' t]);
t_is(full(Gfvv), num_Gfvv, 2, ['Gfvv' t]);

t_is(full(Gtaa), num_Gtaa, 2, ['Gtaa' t]);
t_is(full(Gtav), num_Gtav, 2, ['Gtav' t]);
t_is(full(Gtva), num_Gtva, 2, ['Gtva' t]);
t_is(full(Gtvv), num_Gtvv, 2, ['Gtvv' t]);

%%-----  check d2AIbr_dV2 code  -----
t = ' - d2AIbr_dV2 (squared current magnitudes)';
lam = 10 * rand(nl, 1);
% lam = [1; zeros(nl-1, 1)];
num_Gfaa = zeros(nb, nb);
num_Gfav = zeros(nb, nb);
num_Gfva = zeros(nb, nb);
num_Gfvv = zeros(nb, nb);
num_Gtaa = zeros(nb, nb);
num_Gtav = zeros(nb, nb);
num_Gtva = zeros(nb, nb);
num_Gtvv = zeros(nb, nb);
[dIf_dVa, dIf_dVm, dIt_dVa, dIt_dVm, If, It] = dIbr_dV(branch, Yf, Yt, V);
[dAf_dVa, dAf_dVm, dAt_dVa, dAt_dVm] = ...
                        dAbr_dV(dIf_dVa, dIf_dVm, dIt_dVa, dIt_dVm, If, It);
[Gfaa, Gfav, Gfva, Gfvv] = d2AIbr_dV2(dIf_dVa, dIf_dVm, If, Yf, V, lam);
[Gtaa, Gtav, Gtva, Gtvv] = d2AIbr_dV2(dIt_dVa, dIt_dVm, It, Yt, V, lam);
for i = 1:nb
    Vap = V;
    Vap(i) = Vm(i) * exp(1j * (Va(i) + pert));
    [dIf_dVa_ap, dIf_dVm_ap, dIt_dVa_ap, dIt_dVm_ap, If_ap, It_ap] = ...
        dIbr_dV(branch, Yf, Yt, Vap);
    [dAf_dVa_ap, dAf_dVm_ap, dAt_dVa_ap, dAt_dVm_ap] = ...
        dAbr_dV(dIf_dVa_ap, dIf_dVm_ap, dIt_dVa_ap, dIt_dVm_ap, If_ap, It_ap);
    num_Gfaa(:, i) = (dAf_dVa_ap - dAf_dVa).' * lam / pert;
    num_Gfva(:, i) = (dAf_dVm_ap - dAf_dVm).' * lam / pert;
    num_Gtaa(:, i) = (dAt_dVa_ap - dAt_dVa).' * lam / pert;
    num_Gtva(:, i) = (dAt_dVm_ap - dAt_dVm).' * lam / pert;

    Vmp = V;
    Vmp(i) = (Vm(i) + pert) * exp(1j * Va(i));
    [dIf_dVa_mp, dIf_dVm_mp, dIt_dVa_mp, dIt_dVm_mp, If_mp, It_mp] = ...
        dIbr_dV(branch, Yf, Yt, Vmp);
    [dAf_dVa_mp, dAf_dVm_mp, dAt_dVa_mp, dAt_dVm_mp] = ...
        dAbr_dV(dIf_dVa_mp, dIf_dVm_mp, dIt_dVa_mp, dIt_dVm_mp, If_mp, It_mp);
    num_Gfav(:, i) = (dAf_dVa_mp - dAf_dVa).' * lam / pert;
    num_Gfvv(:, i) = (dAf_dVm_mp - dAf_dVm).' * lam / pert;
    num_Gtav(:, i) = (dAt_dVa_mp - dAt_dVa).' * lam / pert;
    num_Gtvv(:, i) = (dAt_dVm_mp - dAt_dVm).' * lam / pert;
end

t_is(full(Gfaa), num_Gfaa, 3, ['Gfaa' t]);
t_is(full(Gfav), num_Gfav, 3, ['Gfav' t]);
t_is(full(Gfva), num_Gfva, 3, ['Gfva' t]);
t_is(full(Gfvv), num_Gfvv, 2, ['Gfvv' t]);

t_is(full(Gtaa), num_Gtaa, 3, ['Gtaa' t]);
t_is(full(Gtav), num_Gtav, 3, ['Gtav' t]);
t_is(full(Gtva), num_Gtva, 3, ['Gtva' t]);
t_is(full(Gtvv), num_Gtvv, 2, ['Gtvv' t]);

t_end;
