function [Bp, Bpp] = makeB(baseMVA, bus, branch, alg)
%MAKEB   Builds the FDPF matrices, B prime and B double prime.
%   [Bp, Bpp] = makeB(baseMVA, bus, branch, alg) returns the two
%   matrices B prime and B double prime used in the fast decoupled power
%   flow. Does appropriate conversions to p.u.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2003 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%% constants
nb = size(bus, 1);			%% number of buses
nl = size(branch, 1);		%% number of lines

%% define named indices into bus, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
	RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch;

%%-----  form Bp (B prime)  -----
temp_branch = branch;						%% modify a copy of branch
temp_bus = bus;								%% modify a copy of bus
temp_bus(:, BS) = zeros(nb, 1);				%% zero out shunts at buses
temp_branch(:, BR_B) = zeros(nl, 1);		%% zero out line charging shunts
temp_branch(:, TAP) = ones(nl, 1);			%% cancel out taps
if alg == 2									%% if XB method
	temp_branch(:, BR_R) = zeros(nl, 1);		%% zero out line resistance
end
Bp = -imag( makeYbus(baseMVA, temp_bus, temp_branch) );

%%-----  form Bpp (B double prime)  -----
if nargout == 2
	temp_branch = branch;						%% modify a copy of branch
	temp_branch(:, SHIFT) = zeros(nl, 1);		%% zero out phase shifters
	if alg == 3									%% if BX method
		temp_branch(:, BR_R) = zeros(nl, 1);		%% zero out line resistance
	end
	Bpp = -imag( makeYbus(baseMVA, bus, temp_branch) );
end

return;
