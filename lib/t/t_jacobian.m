function t_jacobian(quiet)
%T_JACOBIAN  Numerical tests of partial derivative code.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

if nargin < 1
	quiet = 0;
end

t_begin(12, quiet);

casefile = 'case30';

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
	RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch;

%% run powerflow to get solved case
opt = mpoption('VERBOSE', 0, 'OUT_ALL', 0);
[baseMVA, bus, gen, branch, success, et] = runpf(casefile, opt);

%% switch to internal bus numbering and build admittance matrices
[i2e, bus, gen, branch] = ext2int(bus, gen, branch);
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
Ybus_full	= full(Ybus);
Yf_full 	= full(Yf);
Yt_full 	= full(Yt);
j = sqrt(-1);
V = bus(:, VM) .* exp(j * pi/180 * bus(:, VA));
Vm = abs(V);
Va = angle(V);
f = branch(:, F_BUS);		%% list of "from" buses
t = branch(:, T_BUS);		%% list of "to" buses
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
Vmp = (Vm*ones(1,nb) + pert*eye(nb,nb)) .* (exp(j * Va) * ones(1,nb));
num_dSbus_dVm = full( (Vmp .* conj(Ybus * Vmp) - V*ones(1,nb) .* conj(Ybus * V*ones(1,nb))) / pert );
Vap = (Vm*ones(1,nb)) .* (exp(j * (Va*ones(1,nb) + pert*eye(nb,nb))));
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
Vmp = (Vm*ones(1,nb) + pert*eye(nb,nb)) .* (exp(j * Va) * ones(1,nb));
Vap = (Vm*ones(1,nb)) .* (exp(j * (Va*ones(1,nb) + pert*eye(nb,nb))));
Vmpf = Vmp(f,:);
Vapf = Vap(f,:);
Vmpt = Vmp(t,:);
Vapt = Vap(t,:);
num_dSf_dVm = full( (Vmpf .* conj(Yf * Vmp) - V(f)*ones(1,nb) .* conj(Yf * V*ones(1,nb))) / pert );
num_dSf_dVa = full( (Vapf .* conj(Yf * Vap) - V(f)*ones(1,nb) .* conj(Yf * V*ones(1,nb))) / pert );
num_dSt_dVm = full( (Vmpt .* conj(Yt * Vmp) - V(t)*ones(1,nb) .* conj(Yt * V*ones(1,nb))) / pert );
num_dSt_dVa = full( (Vapt .* conj(Yt * Vap) - V(t)*ones(1,nb) .* conj(Yt * V*ones(1,nb))) / pert );

t_is(dSf_dVm_sp, num_dSf_dVm, 5, 'dSf_dVm (sparse)');
t_is(dSf_dVa_sp, num_dSf_dVa, 5, 'dSf_dVa (sparse)');
t_is(dSt_dVm_sp, num_dSt_dVm, 5, 'dSt_dVm (sparse)');
t_is(dSt_dVa_sp, num_dSt_dVa, 5, 'dSt_dVa (sparse)');
t_is(dSf_dVm_full, num_dSf_dVm, 5, 'dSf_dVm (full)');
t_is(dSf_dVa_full, num_dSf_dVa, 5, 'dSf_dVa (full)');
t_is(dSt_dVm_full, num_dSt_dVm, 5, 'dSt_dVm (full)');
t_is(dSt_dVa_full, num_dSt_dVa, 5, 'dSt_dVa (full)');

t_end;

return;