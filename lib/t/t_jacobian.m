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

t_begin(28, quiet);

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
Ybus_full   = full(Ybus);
Yf_full     = full(Yf);
Yt_full     = full(Yt);
Vm = bus(:, VM);
Va = bus(:, VA) * pi/180;
V = Vm .* exp(1j * Va);
f = branch(:, F_BUS);       %% list of "from" buses
t = branch(:, T_BUS);       %% list of "to" buses
nl = length(f);
nb = length(V);
pert = 1e-8;

%%-----  check dSbus_dV code  -----
%% full matrices
[dSbus_dVm_full, dSbus_dVa_full] = dSbus_dV(Ybus_full, V);

%% sparse matrices
[dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V);
dSbus_dVm_sp = full(dSbus_dVm);
dSbus_dVa_sp = full(dSbus_dVa);

%% compute numerically to compare
Vmp = (Vm*ones(1,nb) + pert*eye(nb,nb)) .* (exp(1j * Va) * ones(1,nb));
Vap = (Vm*ones(1,nb)) .* (exp(1j * (Va*ones(1,nb) + pert*eye(nb,nb))));
num_dSbus_dVm = full( (Vmp .* conj(Ybus * Vmp) - V*ones(1,nb) .* conj(Ybus * V*ones(1,nb))) / pert );
num_dSbus_dVa = full( (Vap .* conj(Ybus * Vap) - V*ones(1,nb) .* conj(Ybus * V*ones(1,nb))) / pert );

t_is(dSbus_dVm_sp, num_dSbus_dVm, 5, 'dSbus_dVm (sparse)');
t_is(dSbus_dVa_sp, num_dSbus_dVa, 5, 'dSbus_dVa (sparse)');
t_is(dSbus_dVm_full, num_dSbus_dVm, 5, 'dSbus_dVm (full)');
t_is(dSbus_dVa_full, num_dSbus_dVa, 5, 'dSbus_dVa (full)');

%%-----  check dSbr_dV code  -----
%% full matrices
[dSf_dVa_full, dSf_dVm_full, dSt_dVa_full, dSt_dVm_full, Sf, St] = dSbr_dV(branch, Yf_full, Yt_full, V);

%% sparse matrices
[dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = dSbr_dV(branch, Yf, Yt, V);
dSf_dVa_sp = full(dSf_dVa);
dSf_dVm_sp = full(dSf_dVm);
dSt_dVa_sp = full(dSt_dVa);
dSt_dVm_sp = full(dSt_dVm);

%% compute numerically to compare
Vmpf = Vmp(f,:);
Vapf = Vap(f,:);
Vmpt = Vmp(t,:);
Vapt = Vap(t,:);
Sf2 = (V(f)*ones(1,nb)) .* conj(Yf * V*ones(1,nb));
St2 = (V(t)*ones(1,nb)) .* conj(Yt * V*ones(1,nb));
Smpf = Vmpf .* conj(Yf * Vmp);
Sapf = Vapf .* conj(Yf * Vap);
Smpt = Vmpt .* conj(Yt * Vmp);
Sapt = Vapt .* conj(Yt * Vap);

num_dSf_dVm = full( (Smpf - Sf2) / pert );
num_dSf_dVa = full( (Sapf - Sf2) / pert );
num_dSt_dVm = full( (Smpt - St2) / pert );
num_dSt_dVa = full( (Sapt - St2) / pert );

t_is(dSf_dVm_sp, num_dSf_dVm, 5, 'dSf_dVm (sparse)');
t_is(dSf_dVa_sp, num_dSf_dVa, 5, 'dSf_dVa (sparse)');
t_is(dSt_dVm_sp, num_dSt_dVm, 5, 'dSt_dVm (sparse)');
t_is(dSt_dVa_sp, num_dSt_dVa, 5, 'dSt_dVa (sparse)');
t_is(dSf_dVm_full, num_dSf_dVm, 5, 'dSf_dVm (full)');
t_is(dSf_dVa_full, num_dSf_dVa, 5, 'dSf_dVa (full)');
t_is(dSt_dVm_full, num_dSt_dVm, 5, 'dSt_dVm (full)');
t_is(dSt_dVa_full, num_dSt_dVa, 5, 'dSt_dVa (full)');

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
num_dAf_dVm = full( (abs(Smpf).^2 - abs(Sf2).^2) / pert );
num_dAf_dVa = full( (abs(Sapf).^2 - abs(Sf2).^2) / pert );
num_dAt_dVm = full( (abs(Smpt).^2 - abs(St2).^2) / pert );
num_dAt_dVa = full( (abs(Sapt).^2 - abs(St2).^2) / pert );

t_is(dAf_dVm_sp, num_dAf_dVm, 4, 'dAf_dVm (sparse)');
t_is(dAf_dVa_sp, num_dAf_dVa, 4, 'dAf_dVa (sparse)');
t_is(dAt_dVm_sp, num_dAt_dVm, 4, 'dAt_dVm (sparse)');
t_is(dAt_dVa_sp, num_dAt_dVa, 4, 'dAt_dVa (sparse)');
t_is(dAf_dVm_full, num_dAf_dVm, 4, 'dAf_dVm (full)');
t_is(dAf_dVa_full, num_dAf_dVa, 4, 'dAf_dVa (full)');
t_is(dAt_dVm_full, num_dAt_dVm, 4, 'dAt_dVm (full)');
t_is(dAt_dVa_full, num_dAt_dVa, 4, 'dAt_dVa (full)');

%%-----  check dIbr_dV code  -----
%% full matrices
[dIf_dVa_full, dIf_dVm_full, dIt_dVa_full, dIt_dVm_full, If, It] = dIbr_dV(branch, Yf_full, Yt_full, V);

%% sparse matrices
[dIf_dVa, dIf_dVm, dIt_dVa, dIt_dVm, If, It] = dIbr_dV(branch, Yf, Yt, V);
dIf_dVa_sp = full(dIf_dVa);
dIf_dVm_sp = full(dIf_dVm);
dIt_dVa_sp = full(dIt_dVa);
dIt_dVm_sp = full(dIt_dVm);

%% compute numerically to compare
num_dIf_dVm = full( (Yf * Vmp - Yf * V*ones(1,nb)) / pert );
num_dIf_dVa = full( (Yf * Vap - Yf * V*ones(1,nb)) / pert );
num_dIt_dVm = full( (Yt * Vmp - Yt * V*ones(1,nb)) / pert );
num_dIt_dVa = full( (Yt * Vap - Yt * V*ones(1,nb)) / pert );

t_is(dIf_dVm_sp, num_dIf_dVm, 5, 'dIf_dVm (sparse)');
t_is(dIf_dVa_sp, num_dIf_dVa, 5, 'dIf_dVa (sparse)');
t_is(dIt_dVm_sp, num_dIt_dVm, 5, 'dIt_dVm (sparse)');
t_is(dIt_dVa_sp, num_dIt_dVa, 5, 'dIt_dVa (sparse)');
t_is(dIf_dVm_full, num_dIf_dVm, 5, 'dIf_dVm (full)');
t_is(dIf_dVa_full, num_dIf_dVa, 5, 'dIf_dVa (full)');
t_is(dIt_dVm_full, num_dIt_dVm, 5, 'dIt_dVm (full)');
t_is(dIt_dVa_full, num_dIt_dVa, 5, 'dIt_dVa (full)');

t_end;
