function [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch)
%MAKEYBUS   Builds the bus admittance matrix and branch admittance matrices.
%   [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch) returns the full
%   bus admittance matrix (i.e. for all buses) and the matrices Yf and Yt
%   which, when multiplied by a complex voltage vector, yield the vector
%   currents injected into each line from the "from" and "to" buses
%   respectively of each line. Does appropriate conversions to p.u.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2003 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%% constants
j = sqrt(-1);
nb = size(bus, 1);			%% number of buses
nl = size(branch, 1);		%% number of lines

%% define named indices into bus, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
	RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch;

%% check that bus numbers are equal to indices to bus (one set of bus numbers)
if any(bus(:, BUS_I) ~= [1:nb]')
	error('buses must appear in order by bus number')
end

%% for each branch, compute the elements of the branch admittance matrix where
%%
%%		| If |   | Yff  Yft |   | Vf |
%%		|    | = |          | * |    |
%%		| It |   | Ytf  Ytt |   | Vt |
%%
stat = branch(:, BR_STATUS);					%% ones at in-service branches
Ys = stat ./ (branch(:, BR_R) + j * branch(:, BR_X));	%% series admittance
Bc = stat .* branch(:, BR_B);							%% line charging susceptance
tap = ones(nl, 1);								%% default tap ratio = 1
i = find(branch(:, TAP));						%% indices of non-zero tap ratios
tap(i) = branch(i, TAP);						%% assign non-zero tap ratios
tap = tap .* exp(-j*pi/180 * branch(:, SHIFT));	%% add phase shifters
Ytt = Ys + j*Bc/2;
Yff = Ytt ./ (tap .* conj(tap));
Yft = - Ys ./ conj(tap);
Ytf = - Ys ./ tap;

%% compute shunt admittance
%% if Psh is the real power consumed by the shunt at V = 1.0 p.u.
%% and Qsh is the reactive power injected by the shunt at V = 1.0 p.u.
%% then Psh - j Qsh = V * conj(Ysh * V) = conj(Ysh) = Gs - j Bs,
%% i.e. Ysh = Psh + j Qsh, so ...
Ysh = (bus(:, GS) + j * bus(:, BS)) / baseMVA;	%% vector of shunt admittances

%% build Ybus
f = branch(:, F_BUS);							%% list of "from" buses
t = branch(:, T_BUS);							%% list of "to" buses
Cf = sparse(f, 1:nl, ones(nl, 1), nb, nl);		%% connection matrix for line & from buses
Ct = sparse(t, 1:nl, ones(nl, 1), nb, nl);		%% connection matrix for line & to buses
Ybus = spdiags(Ysh, 0, nb, nb) + ...			%% shunt admittance
	Cf * spdiags(Yff, 0, nl, nl) * Cf' + ...	%% Yff term of branch admittance
	Cf * spdiags(Yft, 0, nl, nl) * Ct' + ...	%% Yft term of branch admittance
	Ct * spdiags(Ytf, 0, nl, nl) * Cf' + ...	%% Ytf term of branch admittance
	Ct * spdiags(Ytt, 0, nl, nl) * Ct';			%% Ytt term of branch admittance

%% Build Yf and Yt such that Yf * V is the vector of complex branch currents injected
%% at each branch's "from" bus, and Yt is the same for the "to" bus end
if nargout > 1
	i = [[1:nl]'; [1:nl]'];		%% double set of row indices	
	Yf = sparse(i, [f; t], [Yff; Yft]);
	Yt = sparse(i, [f; t], [Ytf; Ytt]);
end

return;
