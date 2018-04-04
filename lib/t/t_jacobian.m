function t_jacobian(quiet)
%T_JACOBIAN  Numerical tests of partial derivative code.

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

%%-----  polar coordinate voltages  -----
vcart = 0;

%%-----  create perturbed voltages  -----
Vap = (Vm*ones(1,nb)) .* (exp(1j * (Va*ones(1,nb) + pert*eye(nb,nb))));
Vmp = (Vm*ones(1,nb) + pert*eye(nb,nb)) .* (exp(1j * Va) * ones(1,nb));

%%-----  check dImis_dV code  -----
%% full matrices
[dImis_dVa_full, dImis_dVm_full] = dImis_dV(Sbus, Ybus_full, V, vcart);

%% sparse matrices
[dImis_dVa, dImis_dVm] = dImis_dV(Sbus, Ybus, V, vcart);
dImis_dVa_sp = full(dImis_dVa);
dImis_dVm_sp = full(dImis_dVm);

%% compute numerically to compare
num_dImis_dVa = full( ( Ybus * Vap - conj(SS./Vap) - Ybus*VV + conj(SS./VV)) / pert );
num_dImis_dVm = full( ( Ybus * Vmp - conj(SS./Vmp) - Ybus*VV + conj(SS./VV)) / pert );

t_is(dImis_dVa_sp, num_dImis_dVa, 5, 'dImis_dVa (sparse)');
t_is(dImis_dVm_sp, num_dImis_dVm, 5, 'dImis_dVm (sparse)');
t_is(dImis_dVa_full, num_dImis_dVa, 5, 'dImis_dVa (full)');
t_is(dImis_dVm_full, num_dImis_dVm, 5, 'dImis_dVm (full)');

%%-----  check dSbus_dV code  -----
%% full matrices
[dSbus_dVm_full, dSbus_dVa_full] = dSbus_dV(Ybus_full, V, vcart);

%% sparse matrices
[dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V, vcart);
dSbus_dVa_sp = full(dSbus_dVa);
dSbus_dVm_sp = full(dSbus_dVm);

%% compute numerically to compare
num_dSbus_dVa = full( (Vap .* conj(Ybus * Vap) - VV .* conj(Ybus * VV)) / pert );
num_dSbus_dVm = full( (Vmp .* conj(Ybus * Vmp) - VV .* conj(Ybus * VV)) / pert );

t_is(dSbus_dVa_sp, num_dSbus_dVa, 5, 'dSbus_dVa (sparse)');
t_is(dSbus_dVm_sp, num_dSbus_dVm, 5, 'dSbus_dVm (sparse)');
t_is(dSbus_dVa_full, num_dSbus_dVa, 5, 'dSbus_dVa (full)');
t_is(dSbus_dVm_full, num_dSbus_dVm, 5, 'dSbus_dVm (full)');

%%-----  check dIbr_dV code  -----
%% full matrices
[dIf_dVa_full, dIf_dVm_full, dIt_dVa_full, dIt_dVm_full, If, It] = dIbr_dV(branch, Yf_full, Yt_full, V, vcart);

%% sparse matrices
[dIf_dVa, dIf_dVm, dIt_dVa, dIt_dVm, If, It] = dIbr_dV(branch, Yf, Yt, V, vcart);
dIf_dVa_sp = full(dIf_dVa);
dIf_dVm_sp = full(dIf_dVm);
dIt_dVa_sp = full(dIt_dVa);
dIt_dVm_sp = full(dIt_dVm);

%% compute numerically to compare
num_dIf_dVa = full( (Yf * Vap - Yf * VV) / pert );
num_dIf_dVm = full( (Yf * Vmp - Yf * VV) / pert );
num_dIt_dVa = full( (Yt * Vap - Yt * VV) / pert );
num_dIt_dVm = full( (Yt * Vmp - Yt * VV) / pert );

t_is(dIf_dVa_sp, num_dIf_dVa, 5, 'dIf_dVa (sparse)');
t_is(dIf_dVm_sp, num_dIf_dVm, 5, 'dIf_dVm (sparse)');
t_is(dIt_dVa_sp, num_dIt_dVa, 5, 'dIt_dVa (sparse)');
t_is(dIt_dVm_sp, num_dIt_dVm, 5, 'dIt_dVm (sparse)');
t_is(dIf_dVa_full, num_dIf_dVa, 5, 'dIf_dVa (full)');
t_is(dIf_dVm_full, num_dIf_dVm, 5, 'dIf_dVm (full)');
t_is(dIt_dVa_full, num_dIt_dVa, 5, 'dIt_dVa (full)');
t_is(dIt_dVm_full, num_dIt_dVm, 5, 'dIt_dVm (full)');

%%-----  check dSbr_dV code  -----
%% full matrices
[dSf_dVa_full, dSf_dVm_full, dSt_dVa_full, dSt_dVm_full, Sf, St] = dSbr_dV(branch, Yf_full, Yt_full, V, vcart);

%% sparse matrices
[dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
dSf_dVa_sp = full(dSf_dVa);
dSf_dVm_sp = full(dSf_dVm);
dSt_dVa_sp = full(dSt_dVa);
dSt_dVm_sp = full(dSt_dVm);

%% compute numerically to compare
Vapf = Vap(f,:);
Vmpf = Vmp(f,:);
Vapt = Vap(t,:);
Vmpt = Vmp(t,:);
Sf2 = (V(f)*ones(1,nb)) .* conj(Yf * VV);
St2 = (V(t)*ones(1,nb)) .* conj(Yt * VV);
Sapf = Vapf .* conj(Yf * Vap);
Smpf = Vmpf .* conj(Yf * Vmp);
Sapt = Vapt .* conj(Yt * Vap);
Smpt = Vmpt .* conj(Yt * Vmp);

num_dSf_dVa = full( (Sapf - Sf2) / pert );
num_dSf_dVm = full( (Smpf - Sf2) / pert );
num_dSt_dVa = full( (Sapt - St2) / pert );
num_dSt_dVm = full( (Smpt - St2) / pert );

t_is(dSf_dVa_sp, num_dSf_dVa, 5, 'dSf_dVa (sparse)');
t_is(dSf_dVm_sp, num_dSf_dVm, 5, 'dSf_dVm (sparse)');
t_is(dSt_dVa_sp, num_dSt_dVa, 5, 'dSt_dVa (sparse)');
t_is(dSt_dVm_sp, num_dSt_dVm, 5, 'dSt_dVm (sparse)');
t_is(dSf_dVa_full, num_dSf_dVa, 5, 'dSf_dVa (full)');
t_is(dSf_dVm_full, num_dSf_dVm, 5, 'dSf_dVm (full)');
t_is(dSt_dVa_full, num_dSt_dVa, 5, 'dSt_dVa (full)');
t_is(dSt_dVm_full, num_dSt_dVm, 5, 'dSt_dVm (full)');

%%-----  check dAbr_dV code  -----
%% full matrices
[dAf_dVa_full, dAf_dVm_full, dAt_dVa_full, dAt_dVm_full] = ...
                        dAbr_dV(dSf_dVa_full, dSf_dVm_full, dSt_dVa_full, dSt_dVm_full, Sf, St);
%% sparse matrices
[dAf_dVa, dAf_dVm, dAt_dVa, dAt_dVm] = ...
                        dAbr_dV(dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St);
dAf_dVa_sp = full(dAf_dVa);
dAf_dVm_sp = full(dAf_dVm);
dAt_dVa_sp = full(dAt_dVa);
dAt_dVm_sp = full(dAt_dVm);

%% compute numerically to compare
num_dAf_dVa = full( (abs(Sapf).^2 - abs(Sf2).^2) / pert );
num_dAf_dVm = full( (abs(Smpf).^2 - abs(Sf2).^2) / pert );
num_dAt_dVa = full( (abs(Sapt).^2 - abs(St2).^2) / pert );
num_dAt_dVm = full( (abs(Smpt).^2 - abs(St2).^2) / pert );

t_is(dAf_dVa_sp, num_dAf_dVa, 4, 'dAf_dVa (sparse)');
t_is(dAf_dVm_sp, num_dAf_dVm, 4, 'dAf_dVm (sparse)');
t_is(dAt_dVa_sp, num_dAt_dVa, 4, 'dAt_dVa (sparse)');
t_is(dAt_dVm_sp, num_dAt_dVm, 4, 'dAt_dVm (sparse)');
t_is(dAf_dVa_full, num_dAf_dVa, 4, 'dAf_dVa (full)');
t_is(dAf_dVm_full, num_dAf_dVm, 4, 'dAf_dVm (full)');
t_is(dAt_dVa_full, num_dAt_dVa, 4, 'dAt_dVa (full)');
t_is(dAt_dVm_full, num_dAt_dVm, 4, 'dAt_dVm (full)');


%%-----  cartesian coordinate voltages  -----
vcart = 1;

%%-----  create perturbed voltages  -----
Vrp = (Vr*ones(1,nb) + pert*eye(nb,nb)) + 1j * Vi * ones(1,nb);
Vip = Vr*ones(1,nb) + 1j * (Vi*ones(1,nb) + pert*eye(nb,nb));

%%-----  check dImis_dV code  -----
%% full matrices
[dImis_dVr_full, dImis_dVi_full] = dImis_dV(Sbus, Ybus_full, V, vcart);

%% sparse matrices
[dImis_dVr, dImis_dVi] = dImis_dV(Sbus, Ybus, V, vcart);
dImis_dVr_sp = full(dImis_dVr);
dImis_dVi_sp = full(dImis_dVi);

%% compute numerically to compare
num_dImis_dVr = full( ( Ybus * Vrp - conj(SS./Vrp) - Ybus*VV + conj(SS./VV)) / pert );
num_dImis_dVi = full( ( Ybus * Vip - conj(SS./Vip) - Ybus*VV + conj(SS./VV)) / pert );

t_is(dImis_dVr_sp, num_dImis_dVr, 5, 'dImis_dVr (sparse)');
t_is(dImis_dVi_sp, num_dImis_dVi, 5, 'dImis_dVi (sparse)');
t_is(dImis_dVr_full, num_dImis_dVr, 5, 'dImis_dVr (full)');
t_is(dImis_dVi_full, num_dImis_dVi, 5, 'dImis_dVi (full)');

%%-----  check dSbus_dV code  -----
%% full matrices
[dSbus_dVr_full, dSbus_dVi_full] = dSbus_dV(Ybus_full, V, vcart);

%% sparse matrices
[dSbus_dVr, dSbus_dVi] = dSbus_dV(Ybus, V, vcart);
dSbus_dVr_sp = full(dSbus_dVr);
dSbus_dVi_sp = full(dSbus_dVi);

%% compute numerically to compare
num_dSbus_dVr = full( (Vrp .* conj(Ybus * Vrp) - VV .* conj(Ybus * VV)) / pert );
num_dSbus_dVi = full( (Vip .* conj(Ybus * Vip) - VV .* conj(Ybus * VV)) / pert );

t_is(dSbus_dVr_sp, num_dSbus_dVr, 5, 'dSbus_dVr (sparse)');
t_is(dSbus_dVi_sp, num_dSbus_dVi, 5, 'dSbus_dVi (sparse)');
t_is(dSbus_dVr_full, num_dSbus_dVr, 5, 'dSbus_dVr (full)');
t_is(dSbus_dVi_full, num_dSbus_dVi, 5, 'dSbus_dVi (full)');

%%-----  check dIbr_dV code  -----
%% full matrices
[dIf_dVr_full, dIf_dVi_full, dIt_dVr_full, dIt_dVi_full, If, It] = dIbr_dV(branch, Yf_full, Yt_full, V, vcart);

%% sparse matrices
[dIf_dVr, dIf_dVi, dIt_dVr, dIt_dVi, If, It] = dIbr_dV(branch, Yf, Yt, V, vcart);
dIf_dVi_sp = full(dIf_dVi);
dIf_dVr_sp = full(dIf_dVr);
dIt_dVi_sp = full(dIt_dVi);
dIt_dVr_sp = full(dIt_dVr);

%% compute numerically to compare
num_dIf_dVr = full( (Yf * Vrp - Yf * VV) / pert );
num_dIf_dVi = full( (Yf * Vip - Yf * VV) / pert );
num_dIt_dVr = full( (Yt * Vrp - Yt * VV) / pert );
num_dIt_dVi = full( (Yt * Vip - Yt * VV) / pert );

t_is(dIf_dVr_sp, num_dIf_dVr, 5, 'dIf_dVr (sparse)');
t_is(dIf_dVi_sp, num_dIf_dVi, 5, 'dIf_dVi (sparse)');
t_is(dIt_dVr_sp, num_dIt_dVr, 5, 'dIt_dVr (sparse)');
t_is(dIt_dVi_sp, num_dIt_dVi, 5, 'dIt_dVi (sparse)');
t_is(dIf_dVr_full, num_dIf_dVr, 5, 'dIf_dVr (full)');
t_is(dIf_dVi_full, num_dIf_dVi, 5, 'dIf_dVi (full)');
t_is(dIt_dVr_full, num_dIt_dVr, 5, 'dIt_dVr (full)');
t_is(dIt_dVi_full, num_dIt_dVi, 5, 'dIt_dVi (full)');

%%-----  check dSbr_dV code  -----
%% full matrices
[dSf_dVr_full, dSf_dVi_full, dSt_dVr_full, dSt_dVi_full, Sf, St] = dSbr_dV(branch, Yf_full, Yt_full, V, vcart);

%% sparse matrices
[dSf_dVr, dSf_dVi, dSt_dVr, dSt_dVi, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
dSf_dVi_sp = full(dSf_dVi);
dSf_dVr_sp = full(dSf_dVr);
dSt_dVi_sp = full(dSt_dVi);
dSt_dVr_sp = full(dSt_dVr);

%% compute numerically to compare
Vrpf = Vrp(f,:);
Vipf = Vip(f,:);
Vrpt = Vrp(t,:);
Vipt = Vip(t,:);
Sf2 = (V(f)*ones(1,nb)) .* conj(Yf * VV);
St2 = (V(t)*ones(1,nb)) .* conj(Yt * VV);
Srpf = Vrpf .* conj(Yf * Vrp);
Sipf = Vipf .* conj(Yf * Vip);
Srpt = Vrpt .* conj(Yt * Vrp);
Sipt = Vipt .* conj(Yt * Vip);

num_dSf_dVr = full( (Srpf - Sf2) / pert );
num_dSf_dVi = full( (Sipf - Sf2) / pert );
num_dSt_dVr = full( (Srpt - St2) / pert );
num_dSt_dVi = full( (Sipt - St2) / pert );

t_is(dSf_dVr_sp, num_dSf_dVr, 5, 'dSf_dVr (sparse)');
t_is(dSf_dVi_sp, num_dSf_dVi, 5, 'dSf_dVi (sparse)');
t_is(dSt_dVr_sp, num_dSt_dVr, 5, 'dSt_dVr (sparse)');
t_is(dSt_dVi_sp, num_dSt_dVi, 5, 'dSt_dVi (sparse)');
t_is(dSf_dVr_full, num_dSf_dVr, 5, 'dSf_dVr (full)');
t_is(dSf_dVi_full, num_dSf_dVi, 5, 'dSf_dVi (full)');
t_is(dSt_dVr_full, num_dSt_dVr, 5, 'dSt_dVr (full)');
t_is(dSt_dVi_full, num_dSt_dVi, 5, 'dSt_dVi (full)');

%%-----  check dAbr_dV code  -----
%% full matrices
[dAf_dVi_full, dAf_dVr_full, dAt_dVi_full, dAt_dVr_full] = ...
                        dAbr_dV(dSf_dVi_full, dSf_dVr_full, dSt_dVi_full, dSt_dVr_full, Sf, St);
%% sparse matrices
[dAf_dVi, dAf_dVr, dAt_dVi, dAt_dVr] = ...
                        dAbr_dV(dSf_dVi, dSf_dVr, dSt_dVi, dSt_dVr, Sf, St);
dAf_dVi_sp = full(dAf_dVi);
dAf_dVr_sp = full(dAf_dVr);
dAt_dVi_sp = full(dAt_dVi);
dAt_dVr_sp = full(dAt_dVr);

%% compute numerically to compare
num_dAf_dVr = full( (abs(Srpf).^2 - abs(Sf2).^2) / pert );
num_dAf_dVi = full( (abs(Sipf).^2 - abs(Sf2).^2) / pert );
num_dAt_dVr = full( (abs(Srpt).^2 - abs(St2).^2) / pert );
num_dAt_dVi = full( (abs(Sipt).^2 - abs(St2).^2) / pert );

t_is(dAf_dVr_sp, num_dAf_dVr, 4, 'dAf_dVr (sparse)');
t_is(dAf_dVi_sp, num_dAf_dVi, 4, 'dAf_dVi (sparse)');
t_is(dAt_dVr_sp, num_dAt_dVr, 4, 'dAt_dVr (sparse)');
t_is(dAt_dVi_sp, num_dAt_dVi, 4, 'dAt_dVi (sparse)');
t_is(dAf_dVr_full, num_dAf_dVr, 4, 'dAf_dVr (full)');
t_is(dAf_dVi_full, num_dAf_dVi, 4, 'dAf_dVi (full)');
t_is(dAt_dVr_full, num_dAt_dVr, 4, 'dAt_dVr (full)');
t_is(dAt_dVi_full, num_dAt_dVi, 4, 'dAt_dVi (full)');

t_end;
