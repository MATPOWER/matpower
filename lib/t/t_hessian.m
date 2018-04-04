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

%%-----  polar coordinate voltages  -----
vcart = 0;

%%-----  check d2Imis_dV2 code  -----
t = ' - d2Imis_dV2 (complex current injections)';
lam = 10 * rand(nb, 1);
num_Haa = zeros(nb, nb);
num_Hav = zeros(nb, nb);
num_Hva = zeros(nb, nb);
num_Hvv = zeros(nb, nb);
[dImis_dVa, dImis_dVm] = dImis_dV(Sbus, Ybus, V, vcart);
[Haa, Hav, Hva, Hvv] = d2Imis_dV2(Sbus, Ybus, V, lam, vcart);
for i = 1:nb
    Vap = V;
    Vap(i) = Vm(i) * exp(1j * (Va(i) + pert));
    [dImis_dVa_ap, dImis_dVm_ap] = dImis_dV(Sbus, Ybus, Vap, vcart);
    num_Haa(:, i) = (dImis_dVa_ap - dImis_dVa).' * lam / pert;
    num_Hva(:, i) = (dImis_dVm_ap - dImis_dVm).' * lam / pert;

    Vmp = V;
    Vmp(i) = (Vm(i) + pert) * exp(1j * Va(i));
    [dImis_dVa_mp, dImis_dVm_mp] = dImis_dV(Sbus, Ybus, Vmp, vcart);
    num_Hav(:, i) = (dImis_dVa_mp - dImis_dVa).' * lam / pert;
    num_Hvv(:, i) = (dImis_dVm_mp - dImis_dVm).' * lam / pert;
end

t_is(full(Haa), num_Haa, 4, ['Haa' t]);
t_is(full(Hav), num_Hav, 4, ['Hav' t]);
t_is(full(Hva), num_Hva, 4, ['Hva' t]);
t_is(full(Hvv), num_Hvv, 4, ['Hvv' t]);

%%-----  check d2Sbus_dV2 code  -----
t = ' - d2Sbus_dV2 (complex power injections)';
lam = 10 * rand(nb, 1);
num_Haa = zeros(nb, nb);
num_Hav = zeros(nb, nb);
num_Hva = zeros(nb, nb);
num_Hvv = zeros(nb, nb);
[dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V, vcart);
[Haa, Hav, Hva, Hvv] = d2Sbus_dV2(Ybus, V, lam, vcart);
for i = 1:nb
    Vap = V;
    Vap(i) = Vm(i) * exp(1j * (Va(i) + pert));
    [dSbus_dVm_ap, dSbus_dVa_ap] = dSbus_dV(Ybus, Vap, vcart);
    num_Haa(:, i) = (dSbus_dVa_ap - dSbus_dVa).' * lam / pert;
    num_Hva(:, i) = (dSbus_dVm_ap - dSbus_dVm).' * lam / pert;

    Vmp = V;
    Vmp(i) = (Vm(i) + pert) * exp(1j * Va(i));
    [dSbus_dVm_mp, dSbus_dVa_mp] = dSbus_dV(Ybus, Vmp, vcart);
    num_Hav(:, i) = (dSbus_dVa_mp - dSbus_dVa).' * lam / pert;
    num_Hvv(:, i) = (dSbus_dVm_mp - dSbus_dVm).' * lam / pert;
end

t_is(full(Haa), num_Haa, 4, ['Haa' t]);
t_is(full(Hav), num_Hav, 4, ['Hav' t]);
t_is(full(Hva), num_Hva, 4, ['Hva' t]);
t_is(full(Hvv), num_Hvv, 4, ['Hvv' t]);

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
[dIf_dVa, dIf_dVm, dIt_dVa, dIt_dVm, If, It] = dIbr_dV(branch, Yf, Yt, V, vcart);
[Gfaa, Gfav, Gfva, Gfvv] = d2Ibr_dV2(Yf, V, lam, vcart);
[Gtaa, Gtav, Gtva, Gtvv] = d2Ibr_dV2(Yt, V, lam, vcart);
for i = 1:nb
    Vap = V;
    Vap(i) = Vm(i) * exp(1j * (Va(i) + pert));
    [dIf_dVa_ap, dIf_dVm_ap, dIt_dVa_ap, dIt_dVm_ap, If_ap, It_ap] = ...
        dIbr_dV(branch, Yf, Yt, Vap, vcart);
    num_Gfaa(:, i) = (dIf_dVa_ap - dIf_dVa).' * lam / pert;
    num_Gfva(:, i) = (dIf_dVm_ap - dIf_dVm).' * lam / pert;
    num_Gtaa(:, i) = (dIt_dVa_ap - dIt_dVa).' * lam / pert;
    num_Gtva(:, i) = (dIt_dVm_ap - dIt_dVm).' * lam / pert;

    Vmp = V;
    Vmp(i) = (Vm(i) + pert) * exp(1j * Va(i));
    [dIf_dVa_mp, dIf_dVm_mp, dIt_dVa_mp, dIt_dVm_mp, If_mp, It_mp] = ...
        dIbr_dV(branch, Yf, Yt, Vmp, vcart);
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
[dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
[Gfaa, Gfav, Gfva, Gfvv] = d2Sbr_dV2(Cf, Yf, V, lam, vcart);
[Gtaa, Gtav, Gtva, Gtvv] = d2Sbr_dV2(Ct, Yt, V, lam, vcart);
for i = 1:nb
    Vap = V;
    Vap(i) = Vm(i) * exp(1j * (Va(i) + pert));
    [dSf_dVa_ap, dSf_dVm_ap, dSt_dVa_ap, dSt_dVm_ap, Sf_ap, St_ap] = ...
        dSbr_dV(branch, Yf, Yt, Vap, vcart);
    num_Gfaa(:, i) = (dSf_dVa_ap - dSf_dVa).' * lam / pert;
    num_Gfva(:, i) = (dSf_dVm_ap - dSf_dVm).' * lam / pert;
    num_Gtaa(:, i) = (dSt_dVa_ap - dSt_dVa).' * lam / pert;
    num_Gtva(:, i) = (dSt_dVm_ap - dSt_dVm).' * lam / pert;

    Vmp = V;
    Vmp(i) = (Vm(i) + pert) * exp(1j * Va(i));
    [dSf_dVa_mp, dSf_dVm_mp, dSt_dVa_mp, dSt_dVm_mp, Sf_mp, St_mp] = ...
        dSbr_dV(branch, Yf, Yt, Vmp, vcart);
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

%%-----  check d2Abr_dV2 code  -----
t = ' - d2Abr_dV2 (squared current magnitudes)';
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
d2If_dV2 = @(V, mu)d2Ibr_dV2(Yf, V, mu, vcart);
d2It_dV2 = @(V, mu)d2Ibr_dV2(Yt, V, mu, vcart);
[dIf_dVa, dIf_dVm, dIt_dVa, dIt_dVm, If, It] = dIbr_dV(branch, Yf, Yt, V, vcart);
[dAf_dVa, dAf_dVm, dAt_dVa, dAt_dVm] = ...
                        dAbr_dV(dIf_dVa, dIf_dVm, dIt_dVa, dIt_dVm, If, It);
[Gfaa, Gfav, Gfva, Gfvv] = d2Abr_dV2(d2If_dV2, dIf_dVa, dIf_dVm, If, V, lam);
[Gtaa, Gtav, Gtva, Gtvv] = d2Abr_dV2(d2It_dV2, dIt_dVa, dIt_dVm, It, V, lam);
for i = 1:nb
    Vap = V;
    Vap(i) = Vm(i) * exp(1j * (Va(i) + pert));
    [dIf_dVa_ap, dIf_dVm_ap, dIt_dVa_ap, dIt_dVm_ap, If_ap, It_ap] = ...
        dIbr_dV(branch, Yf, Yt, Vap, vcart);
    [dAf_dVa_ap, dAf_dVm_ap, dAt_dVa_ap, dAt_dVm_ap] = ...
        dAbr_dV(dIf_dVa_ap, dIf_dVm_ap, dIt_dVa_ap, dIt_dVm_ap, If_ap, It_ap);
    num_Gfaa(:, i) = (dAf_dVa_ap - dAf_dVa).' * lam / pert;
    num_Gfva(:, i) = (dAf_dVm_ap - dAf_dVm).' * lam / pert;
    num_Gtaa(:, i) = (dAt_dVa_ap - dAt_dVa).' * lam / pert;
    num_Gtva(:, i) = (dAt_dVm_ap - dAt_dVm).' * lam / pert;

    Vmp = V;
    Vmp(i) = (Vm(i) + pert) * exp(1j * Va(i));
    [dIf_dVa_mp, dIf_dVm_mp, dIt_dVa_mp, dIt_dVm_mp, If_mp, It_mp] = ...
        dIbr_dV(branch, Yf, Yt, Vmp, vcart);
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

%%-----  check d2Abr_dV2 code  -----
t = ' - d2Abr_dV2 (squared apparent power flows)';
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
d2Sf_dV2 = @(V, mu)d2Sbr_dV2(Cf, Yf, V, mu, vcart);
d2St_dV2 = @(V, mu)d2Sbr_dV2(Ct, Yt, V, mu, vcart);
[dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
[dAf_dVa, dAf_dVm, dAt_dVa, dAt_dVm] = ...
                        dAbr_dV(dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St);
[Gfaa, Gfav, Gfva, Gfvv] = d2Abr_dV2(d2Sf_dV2, dSf_dVa, dSf_dVm, Sf, V, lam);
[Gtaa, Gtav, Gtva, Gtvv] = d2Abr_dV2(d2St_dV2, dSt_dVa, dSt_dVm, St, V, lam);
for i = 1:nb
    Vap = V;
    Vap(i) = Vm(i) * exp(1j * (Va(i) + pert));
    [dSf_dVa_ap, dSf_dVm_ap, dSt_dVa_ap, dSt_dVm_ap, Sf_ap, St_ap] = ...
        dSbr_dV(branch, Yf, Yt, Vap, vcart);
    [dAf_dVa_ap, dAf_dVm_ap, dAt_dVa_ap, dAt_dVm_ap] = ...
        dAbr_dV(dSf_dVa_ap, dSf_dVm_ap, dSt_dVa_ap, dSt_dVm_ap, Sf_ap, St_ap);
    num_Gfaa(:, i) = (dAf_dVa_ap - dAf_dVa).' * lam / pert;
    num_Gfva(:, i) = (dAf_dVm_ap - dAf_dVm).' * lam / pert;
    num_Gtaa(:, i) = (dAt_dVa_ap - dAt_dVa).' * lam / pert;
    num_Gtva(:, i) = (dAt_dVm_ap - dAt_dVm).' * lam / pert;

    Vmp = V;
    Vmp(i) = (Vm(i) + pert) * exp(1j * Va(i));
    [dSf_dVa_mp, dSf_dVm_mp, dSt_dVa_mp, dSt_dVm_mp, Sf_mp, St_mp] = ...
        dSbr_dV(branch, Yf, Yt, Vmp, vcart);
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

%%-----  check d2Abr_dV2 code  -----
t = ' - d2Abr_dV2 (squared real power flows)';
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
d2Sf_dV2 = @(V, mu)d2Sbr_dV2(Cf, Yf, V, mu, vcart);
d2St_dV2 = @(V, mu)d2Sbr_dV2(Ct, Yt, V, mu, vcart);
[dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
[dAf_dVa, dAf_dVm, dAt_dVa, dAt_dVm] = ...
                        dAbr_dV(real(dSf_dVa), real(dSf_dVm), real(dSt_dVa), real(dSt_dVm), real(Sf), real(St));
[Gfaa, Gfav, Gfva, Gfvv] = d2Abr_dV2(d2Sf_dV2, real(dSf_dVa), real(dSf_dVm), real(Sf), V, lam);
[Gtaa, Gtav, Gtva, Gtvv] = d2Abr_dV2(d2St_dV2, real(dSt_dVa), real(dSt_dVm), real(St), V, lam);
for i = 1:nb
    Vap = V;
    Vap(i) = Vm(i) * exp(1j * (Va(i) + pert));
    [dSf_dVa_ap, dSf_dVm_ap, dSt_dVa_ap, dSt_dVm_ap, Sf_ap, St_ap] = ...
        dSbr_dV(branch, Yf, Yt, Vap, vcart);
    [dAf_dVa_ap, dAf_dVm_ap, dAt_dVa_ap, dAt_dVm_ap] = ...
        dAbr_dV(real(dSf_dVa_ap), real(dSf_dVm_ap), real(dSt_dVa_ap), real(dSt_dVm_ap), real(Sf_ap), real(St_ap));
    num_Gfaa(:, i) = (dAf_dVa_ap - dAf_dVa).' * lam / pert;
    num_Gfva(:, i) = (dAf_dVm_ap - dAf_dVm).' * lam / pert;
    num_Gtaa(:, i) = (dAt_dVa_ap - dAt_dVa).' * lam / pert;
    num_Gtva(:, i) = (dAt_dVm_ap - dAt_dVm).' * lam / pert;

    Vmp = V;
    Vmp(i) = (Vm(i) + pert) * exp(1j * Va(i));
    [dSf_dVa_mp, dSf_dVm_mp, dSt_dVa_mp, dSt_dVm_mp, Sf_mp, St_mp] = ...
        dSbr_dV(branch, Yf, Yt, Vmp, vcart);
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

%%-----  cartesian coordinate voltages  -----
vcart = 1;

%%-----  check d2Imis_dV2 code  -----
t = ' - d2Imis_dV2 (complex current injections)';
lam = 10 * rand(nb, 1);
num_Hrr = zeros(nb, nb);
num_Hri = zeros(nb, nb);
num_Hir = zeros(nb, nb);
num_Hii = zeros(nb, nb);
[dImis_dVr, dImis_dVi] = dImis_dV(Sbus, Ybus, V, vcart);
[Hrr, Hri, Hir, Hii] = d2Imis_dV2(Sbus, Ybus, V, lam, vcart);
for i = 1:nb
    Vrp = V;
    Vrp(i) = (Vr(i) + pert) + 1j * Vi(i);
    [dImis_dVr_rp, dImis_dVi_rp] = dImis_dV(Sbus, Ybus, Vrp, vcart);
    num_Hrr(:, i) = (dImis_dVr_rp - dImis_dVr).' * lam / pert;
    num_Hir(:, i) = (dImis_dVi_rp - dImis_dVi).' * lam / pert;

    Vip = V;
    Vip(i) = Vr(i) + 1j * (Vi(i) + pert);
    [dImis_dVr_ip, dImis_dVi_ip] = dImis_dV(Sbus, Ybus, Vip, vcart);
    num_Hri(:, i) = (dImis_dVr_ip - dImis_dVr).' * lam / pert;
    num_Hii(:, i) = (dImis_dVi_ip - dImis_dVi).' * lam / pert;
end

t_is(full(Hrr), num_Hrr, 4, ['Hrr' t]);
t_is(full(Hri), num_Hri, 4, ['Hri' t]);
t_is(full(Hir), num_Hir, 4, ['Hir' t]);
t_is(full(Hii), num_Hii, 4, ['Hii' t]);

%%-----  check d2Sbus_dV2 code  -----
t = ' - d2Sbus_dV2 (complex power injections)';
lam = 10 * rand(nb, 1);
num_Hrr = zeros(nb, nb);
num_Hri = zeros(nb, nb);
num_Hir = zeros(nb, nb);
num_Hii = zeros(nb, nb);
[dSbus_dVr, dSbus_dVi] = dSbus_dV(Ybus, V, vcart);
[Hrr, Hri, Hir, Hii] = d2Sbus_dV2(Ybus, V, lam, vcart);
for i = 1:nb
    Vrp = V;
    Vrp(i) = (Vr(i) + pert) + 1j * Vi(i);
    [dSbus_dVr_rp, dSbus_dVi_rp] = dSbus_dV(Ybus, Vrp, vcart);
    num_Hrr(:, i) = (dSbus_dVr_rp - dSbus_dVr).' * lam / pert;
    num_Hir(:, i) = (dSbus_dVi_rp - dSbus_dVi).' * lam / pert;

    Vip = V;
    Vip(i) = Vr(i) + 1j * (Vi(i) + pert);
    [dSbus_dVr_ip, dSbus_dVi_ip] = dSbus_dV(Ybus, Vip, vcart);
    num_Hri(:, i) = (dSbus_dVr_ip - dSbus_dVr).' * lam / pert;
    num_Hii(:, i) = (dSbus_dVi_ip - dSbus_dVi).' * lam / pert;
end

t_is(full(Hrr), num_Hrr, 4, ['Hrr' t]);
t_is(full(Hri), num_Hri, 4, ['Hri' t]);
t_is(full(Hir), num_Hir, 4, ['Hir' t]);
t_is(full(Hii), num_Hii, 4, ['Hii' t]);

%%-----  check d2Ibr_dV2 code  -----
t = ' - d2Ibr_dV2 (complex currents)';
lam = 10 * rand(nl, 1);
% lam = [1; zeros(nl-1, 1)];
num_Gfrr = zeros(nb, nb);
num_Gfri = zeros(nb, nb);
num_Gfir = zeros(nb, nb);
num_Gfii = zeros(nb, nb);
num_Gtrr = zeros(nb, nb);
num_Gtri = zeros(nb, nb);
num_Gtir = zeros(nb, nb);
num_Gtii = zeros(nb, nb);
[dIf_dVr, dIf_dVi, dIt_dVr, dIt_dVi, If, It] = dIbr_dV(branch, Yf, Yt, V, vcart);
[Gfrr, Gfri, Gfir, Gfii] = d2Ibr_dV2(Yf, V, lam, vcart);
[Gtrr, Gtri, Gtir, Gtii] = d2Ibr_dV2(Yt, V, lam, vcart);
for i = 1:nb
    Vrp = V;
    Vrp(i) = (Vr(i) + pert) + 1j * Vi(i);
    [dIf_dVr_rp, dIf_dVi_rp, dIt_dVr_rp, dIt_dVi_rp, If_rp, It_rp] = ...
        dIbr_dV(branch, Yf, Yt, Vrp, vcart);
    num_Gfrr(:, i) = (dIf_dVr_rp - dIf_dVr).' * lam / pert;
    num_Gfir(:, i) = (dIf_dVi_rp - dIf_dVi).' * lam / pert;
    num_Gtrr(:, i) = (dIt_dVr_rp - dIt_dVr).' * lam / pert;
    num_Gtir(:, i) = (dIt_dVi_rp - dIt_dVi).' * lam / pert;

    Vip = V;
    Vip(i) = Vr(i) + 1j * (Vi(i) + pert);
    [dIf_dVr_ip, dIf_dVi_ip, dIt_dVr_ip, dIt_dVi_ip, If_ip, It_ip] = ...
        dIbr_dV(branch, Yf, Yt, Vip, vcart);
    num_Gfri(:, i) = (dIf_dVr_ip - dIf_dVr).' * lam / pert;
    num_Gfii(:, i) = (dIf_dVi_ip - dIf_dVi).' * lam / pert;
    num_Gtri(:, i) = (dIt_dVr_ip - dIt_dVr).' * lam / pert;
    num_Gtii(:, i) = (dIt_dVi_ip - dIt_dVi).' * lam / pert;
end

t_is(full(Gfrr), num_Gfrr, 4, ['Gfrr' t]);
t_is(full(Gfri), num_Gfri, 4, ['Gfri' t]);
t_is(full(Gfir), num_Gfir, 4, ['Gfir' t]);
t_is(full(Gfii), num_Gfii, 4, ['Gfii' t]);

t_is(full(Gtrr), num_Gtrr, 4, ['Gtrr' t]);
t_is(full(Gtri), num_Gtri, 4, ['Gtri' t]);
t_is(full(Gtir), num_Gtir, 4, ['Gtir' t]);
t_is(full(Gtii), num_Gtii, 4, ['Gtii' t]);

%%-----  check d2Sbr_dV2 code  -----
t = ' - d2Sbr_dV2 (complex power flows)';
lam = 10 * rand(nl, 1);
% lam = [1; zeros(nl-1, 1)];
num_Gfrr = zeros(nb, nb);
num_Gfri = zeros(nb, nb);
num_Gfir = zeros(nb, nb);
num_Gfii = zeros(nb, nb);
num_Gtrr = zeros(nb, nb);
num_Gtri = zeros(nb, nb);
num_Gtir = zeros(nb, nb);
num_Gtii = zeros(nb, nb);
[dSf_dVr, dSf_dVi, dSt_dVr, dSt_dVi, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
[Gfrr, Gfri, Gfir, Gfii] = d2Sbr_dV2(Cf, Yf, V, lam, vcart);
[Gtrr, Gtri, Gtir, Gtii] = d2Sbr_dV2(Ct, Yt, V, lam, vcart);
for i = 1:nb
    Vrp = V;
    Vrp(i) = (Vr(i) + pert) + 1j * Vi(i);
    [dSf_dVr_rp, dSf_dVi_rp, dSt_dVr_rp, dSt_dVi_rp, Sf_rp, St_rp] = ...
        dSbr_dV(branch, Yf, Yt, Vrp, vcart);
    num_Gfrr(:, i) = (dSf_dVr_rp - dSf_dVr).' * lam / pert;
    num_Gfir(:, i) = (dSf_dVi_rp - dSf_dVi).' * lam / pert;
    num_Gtrr(:, i) = (dSt_dVr_rp - dSt_dVr).' * lam / pert;
    num_Gtir(:, i) = (dSt_dVi_rp - dSt_dVi).' * lam / pert;

    Vip = V;
    Vip(i) = Vr(i) + 1j * (Vi(i) + pert);
    [dSf_dVr_ip, dSf_dVi_ip, dSt_dVr_ip, dSt_dVi_ip, Sf_ip, St_ip] = ...
        dSbr_dV(branch, Yf, Yt, Vip, vcart);
    num_Gfri(:, i) = (dSf_dVr_ip - dSf_dVr).' * lam / pert;
    num_Gfii(:, i) = (dSf_dVi_ip - dSf_dVi).' * lam / pert;
    num_Gtri(:, i) = (dSt_dVr_ip - dSt_dVr).' * lam / pert;
    num_Gtii(:, i) = (dSt_dVi_ip - dSt_dVi).' * lam / pert;
end

t_is(full(Gfrr), num_Gfrr, 4, ['Gfrr' t]);
t_is(full(Gfri), num_Gfri, 4, ['Gfri' t]);
t_is(full(Gfir), num_Gfir, 4, ['Gfir' t]);
t_is(full(Gfii), num_Gfii, 4, ['Gfii' t]);

t_is(full(Gtrr), num_Gtrr, 4, ['Gtrr' t]);
t_is(full(Gtri), num_Gtri, 4, ['Gtri' t]);
t_is(full(Gtir), num_Gtir, 4, ['Gtir' t]);
t_is(full(Gtii), num_Gtii, 4, ['Gtii' t]);

%%-----  check d2Abr_dV2 code  -----
t = ' - d2Abr_dV2 (squared current magnitudes)';
lam = 10 * rand(nl, 1);
% lam = [1; zeros(nl-1, 1)];
num_Gfrr = zeros(nb, nb);
num_Gfri = zeros(nb, nb);
num_Gfir = zeros(nb, nb);
num_Gfii = zeros(nb, nb);
num_Gtrr = zeros(nb, nb);
num_Gtri = zeros(nb, nb);
num_Gtir = zeros(nb, nb);
num_Gtii = zeros(nb, nb);
d2If_dV2 = @(V, mu)d2Ibr_dV2(Yf, V, mu, vcart);
d2It_dV2 = @(V, mu)d2Ibr_dV2(Yt, V, mu, vcart);
[dIf_dVr, dIf_dVi, dIt_dVr, dIt_dVi, If, It] = dIbr_dV(branch, Yf, Yt, V, vcart);
[dAf_dVr, dAf_dVi, dAt_dVr, dAt_dVi] = ...
                        dAbr_dV(dIf_dVr, dIf_dVi, dIt_dVr, dIt_dVi, If, It);
[Gfrr, Gfri, Gfir, Gfii] = d2Abr_dV2(d2If_dV2, dIf_dVr, dIf_dVi, If, V, lam);
[Gtrr, Gtri, Gtir, Gtii] = d2Abr_dV2(d2It_dV2, dIt_dVr, dIt_dVi, It, V, lam);
for i = 1:nb
    Vrp = V;
    Vrp(i) = (Vr(i) + pert) + 1j * Vi(i);
    [dIf_dVr_rp, dIf_dVi_rp, dIt_dVr_rp, dIt_dVi_rp, If_rp, It_rp] = ...
        dIbr_dV(branch, Yf, Yt, Vrp, vcart);
    [dAf_dVr_rp, dAf_dVi_rp, dAt_dVr_rp, dAt_dVi_rp] = ...
        dAbr_dV(dIf_dVr_rp, dIf_dVi_rp, dIt_dVr_rp, dIt_dVi_rp, If_rp, It_rp);
    num_Gfrr(:, i) = (dAf_dVr_rp - dAf_dVr).' * lam / pert;
    num_Gfir(:, i) = (dAf_dVi_rp - dAf_dVi).' * lam / pert;
    num_Gtrr(:, i) = (dAt_dVr_rp - dAt_dVr).' * lam / pert;
    num_Gtir(:, i) = (dAt_dVi_rp - dAt_dVi).' * lam / pert;

    Vip = V;
    Vip(i) = Vr(i) + 1j * (Vi(i) + pert);
    [dIf_dVr_ip, dIf_dVi_ip, dIt_dVr_ip, dIt_dVi_ip, If_ip, It_ip] = ...
        dIbr_dV(branch, Yf, Yt, Vip, vcart);
    [dAf_dVr_ip, dAf_dVi_ip, dAt_dVr_ip, dAt_dVi_ip] = ...
        dAbr_dV(dIf_dVr_ip, dIf_dVi_ip, dIt_dVr_ip, dIt_dVi_ip, If_ip, It_ip);
    num_Gfri(:, i) = (dAf_dVr_ip - dAf_dVr).' * lam / pert;
    num_Gfii(:, i) = (dAf_dVi_ip - dAf_dVi).' * lam / pert;
    num_Gtri(:, i) = (dAt_dVr_ip - dAt_dVr).' * lam / pert;
    num_Gtii(:, i) = (dAt_dVi_ip - dAt_dVi).' * lam / pert;
end

t_is(full(Gfrr), num_Gfrr, 2, ['Gfrr' t]);
t_is(full(Gfri), num_Gfri, 3, ['Gfri' t]);
t_is(full(Gfir), num_Gfir, 3, ['Gfir' t]);
t_is(full(Gfii), num_Gfii, 3, ['Gfii' t]);

t_is(full(Gtrr), num_Gtrr, 2, ['Gtrr' t]);
t_is(full(Gtri), num_Gtri, 3, ['Gtri' t]);
t_is(full(Gtir), num_Gtir, 3, ['Gtir' t]);
t_is(full(Gtii), num_Gtii, 3, ['Gtii' t]);

%%-----  check d2Abr_dV2 code  -----
t = ' - d2Abr_dV2 (squared apparent power flows)';
lam = 10 * rand(nl, 1);
% lam = [1; zeros(nl-1, 1)];
num_Gfrr = zeros(nb, nb);
num_Gfri = zeros(nb, nb);
num_Gfir = zeros(nb, nb);
num_Gfii = zeros(nb, nb);
num_Gtrr = zeros(nb, nb);
num_Gtri = zeros(nb, nb);
num_Gtir = zeros(nb, nb);
num_Gtii = zeros(nb, nb);
d2Sf_dV2 = @(V, mu)d2Sbr_dV2(Cf, Yf, V, mu, vcart);
d2St_dV2 = @(V, mu)d2Sbr_dV2(Ct, Yt, V, mu, vcart);
[dSf_dVr, dSf_dVi, dSt_dVr, dSt_dVi, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
[dAf_dVr, dAf_dVi, dAt_dVr, dAt_dVi] = ...
                        dAbr_dV(dSf_dVr, dSf_dVi, dSt_dVr, dSt_dVi, Sf, St);
[Gfrr, Gfri, Gfir, Gfii] = d2Abr_dV2(d2Sf_dV2, dSf_dVr, dSf_dVi, Sf, V, lam);
[Gtrr, Gtri, Gtir, Gtii] = d2Abr_dV2(d2St_dV2, dSt_dVr, dSt_dVi, St, V, lam);
for i = 1:nb
    Vrp = V;
    Vrp(i) = (Vr(i) + pert) + 1j * Vi(i);
    [dSf_dVr_rp, dSf_dVi_rp, dSt_dVr_rp, dSt_dVi_rp, Sf_rp, St_rp] = ...
        dSbr_dV(branch, Yf, Yt, Vrp, vcart);
    [dAf_dVr_rp, dAf_dVi_rp, dAt_dVr_rp, dAt_dVi_rp] = ...
        dAbr_dV(dSf_dVr_rp, dSf_dVi_rp, dSt_dVr_rp, dSt_dVi_rp, Sf_rp, St_rp);
    num_Gfrr(:, i) = (dAf_dVr_rp - dAf_dVr).' * lam / pert;
    num_Gfir(:, i) = (dAf_dVi_rp - dAf_dVi).' * lam / pert;
    num_Gtrr(:, i) = (dAt_dVr_rp - dAt_dVr).' * lam / pert;
    num_Gtir(:, i) = (dAt_dVi_rp - dAt_dVi).' * lam / pert;

    Vip = V;
    Vip(i) = Vr(i) + 1j * (Vi(i) + pert);
    [dSf_dVr_ip, dSf_dVi_ip, dSt_dVr_ip, dSt_dVi_ip, Sf_ip, St_ip] = ...
        dSbr_dV(branch, Yf, Yt, Vip, vcart);
    [dAf_dVr_ip, dAf_dVi_ip, dAt_dVr_ip, dAt_dVi_ip] = ...
        dAbr_dV(dSf_dVr_ip, dSf_dVi_ip, dSt_dVr_ip, dSt_dVi_ip, Sf_ip, St_ip);
    num_Gfri(:, i) = (dAf_dVr_ip - dAf_dVr).' * lam / pert;
    num_Gfii(:, i) = (dAf_dVi_ip - dAf_dVi).' * lam / pert;
    num_Gtri(:, i) = (dAt_dVr_ip - dAt_dVr).' * lam / pert;
    num_Gtii(:, i) = (dAt_dVi_ip - dAt_dVi).' * lam / pert;
end

t_is(full(Gfrr), num_Gfrr, 2, ['Gfrr' t]);
t_is(full(Gfri), num_Gfri, 2, ['Gfri' t]);
t_is(full(Gfir), num_Gfir, 2, ['Gfir' t]);
t_is(full(Gfii), num_Gfii, 2, ['Gfii' t]);

t_is(full(Gtrr), num_Gtrr, 2, ['Gtrr' t]);
t_is(full(Gtri), num_Gtri, 2, ['Gtri' t]);
t_is(full(Gtir), num_Gtir, 2, ['Gtir' t]);
t_is(full(Gtii), num_Gtii, 2, ['Gtii' t]);

%%-----  check d2Abr_dV2 code  -----
t = ' - d2Abr_dV2 (squared real power flows)';
lam = 10 * rand(nl, 1);
% lam = [1; zeros(nl-1, 1)];
num_Gfrr = zeros(nb, nb);
num_Gfri = zeros(nb, nb);
num_Gfir = zeros(nb, nb);
num_Gfii = zeros(nb, nb);
num_Gtrr = zeros(nb, nb);
num_Gtri = zeros(nb, nb);
num_Gtir = zeros(nb, nb);
num_Gtii = zeros(nb, nb);
d2Sf_dV2 = @(V, mu)d2Sbr_dV2(Cf, Yf, V, mu, vcart);
d2St_dV2 = @(V, mu)d2Sbr_dV2(Ct, Yt, V, mu, vcart);
[dSf_dVr, dSf_dVi, dSt_dVr, dSt_dVi, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
[dAf_dVr, dAf_dVi, dAt_dVr, dAt_dVi] = ...
                        dAbr_dV(real(dSf_dVr), real(dSf_dVi), real(dSt_dVr), real(dSt_dVi), real(Sf), real(St));
[Gfrr, Gfri, Gfir, Gfii] = d2Abr_dV2(d2Sf_dV2, real(dSf_dVr), real(dSf_dVi), real(Sf), V, lam);
[Gtrr, Gtri, Gtir, Gtii] = d2Abr_dV2(d2St_dV2, real(dSt_dVr), real(dSt_dVi), real(St), V, lam);
for i = 1:nb
    Vrp = V;
    Vrp(i) = (Vr(i) + pert) + 1j * Vi(i);
    [dSf_dVr_rp, dSf_dVi_rp, dSt_dVr_rp, dSt_dVi_rp, Sf_rp, St_rp] = ...
        dSbr_dV(branch, Yf, Yt, Vrp, vcart);
    [dAf_dVr_rp, dAf_dVi_rp, dAt_dVr_rp, dAt_dVi_rp] = ...
        dAbr_dV(real(dSf_dVr_rp), real(dSf_dVi_rp), real(dSt_dVr_rp), real(dSt_dVi_rp), real(Sf_rp), real(St_rp));
    num_Gfrr(:, i) = (dAf_dVr_rp - dAf_dVr).' * lam / pert;
    num_Gfir(:, i) = (dAf_dVi_rp - dAf_dVi).' * lam / pert;
    num_Gtrr(:, i) = (dAt_dVr_rp - dAt_dVr).' * lam / pert;
    num_Gtir(:, i) = (dAt_dVi_rp - dAt_dVi).' * lam / pert;

    Vip = V;
    Vip(i) = Vr(i) + 1j * (Vi(i) + pert);
    [dSf_dVr_ip, dSf_dVi_ip, dSt_dVr_ip, dSt_dVi_ip, Sf_ip, St_ip] = ...
        dSbr_dV(branch, Yf, Yt, Vip, vcart);
    [dAf_dVr_ip, dAf_dVi_ip, dAt_dVr_ip, dAt_dVi_ip] = ...
        dAbr_dV(real(dSf_dVr_ip), real(dSf_dVi_ip), real(dSt_dVr_ip), real(dSt_dVi_ip), real(Sf_ip), real(St_ip));
    num_Gfri(:, i) = (dAf_dVr_ip - dAf_dVr).' * lam / pert;
    num_Gfii(:, i) = (dAf_dVi_ip - dAf_dVi).' * lam / pert;
    num_Gtri(:, i) = (dAt_dVr_ip - dAt_dVr).' * lam / pert;
    num_Gtii(:, i) = (dAt_dVi_ip - dAt_dVi).' * lam / pert;
end

t_is(full(Gfrr), num_Gfrr, 2, ['Gfrr' t]);
t_is(full(Gfri), num_Gfri, 2, ['Gfri' t]);
t_is(full(Gfir), num_Gfir, 2, ['Gfir' t]);
t_is(full(Gfii), num_Gfii, 2, ['Gfii' t]);

t_is(full(Gtrr), num_Gtrr, 2, ['Gtrr' t]);
t_is(full(Gtri), num_Gtri, 2, ['Gtri' t]);
t_is(full(Gtir), num_Gtir, 2, ['Gtir' t]);
t_is(full(Gtii), num_Gtii, 2, ['Gtii' t]);

t_end;
