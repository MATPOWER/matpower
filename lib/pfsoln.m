function [bus, gen, branch] = pfsoln(baseMVA, bus0, gen0, branch0, Ybus, Yf, Yt, V, ref, pv, pq);
%PFSOLN  Updates bus, gen, branch data structures to match power flow soln.
%   [bus, gen, branch] = pfsoln(baseMVA, bus0, gen0, branch0, ...
%                                   Ybus, Yf, Yt, V, ref, pv, pq)

%   MATPOWER Version 2.0
%   by Ray Zimmerman, PSERC Cornell    12/19/97
%   Copyright (c) 1996, 1997 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%% constants
j = sqrt(-1);
nl = size(branch0, 1);		%% number of lines

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
	GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
	RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch;

%% initialize return values
bus		= bus0;
gen		= gen0;
branch	= branch0;

%%----- update bus voltages -----
bus(:, VM) = abs(V);
bus(:, VA) = angle(V) * 180 / pi;

%%----- update Qg for all gens and Pg for swing bus -----
%% generator info
on = find(gen(:, GEN_STATUS));				%% which generators are on?
gbus = gen(on, GEN_BUS);					%% what buses are they at?
refgen = find(gen(on, GEN_BUS) == ref);		%% which is the reference gen?

%% compute total injected bus powers
%% This is slow in Matlab 5 ...
% Sg = V(gbus) .* conj(Ybus(gbus, :) * V);
%% ... so we do this instead ...
temp = Ybus.';
Sg = V(gbus) .* conj(temp(:, gbus).' * V);

%% update Qg for all generators
gen(:, QG) = zeros(size(gen, 1), 1);				%% zero out all Qg
gen(on, QG) = imag(Sg) * baseMVA + bus(gbus, QD);	%% inj Q + local Qd

%% update Pg for swing bus
gen(on(refgen), PG) = real(Sg(refgen)) * baseMVA + bus(ref, PD);	%% inj P + local Pd

%%----- update/compute branch power flows -----
Sf = V(branch(:, F_BUS)) .* conj(Yf * V) * baseMVA;	%% complex power at "from" bus
St = V(branch(:, T_BUS)) .* conj(Yt * V) * baseMVA;	%% complex power injected at "to" bus
branch(:, [PF, QF, PT, QT]) = [real(Sf) imag(Sf) real(St) imag(St)];

return;
