function [i2e, bus, gen, branch, area] = ext2int(bus, gen, branch, area)
%EXT2INT   Converts external to internal bus numbering.
%   [i2e, bus, gen, branch, area] = ext2int(bus, gen, branch, area) converts
%   external bus numbers (possibly non-consecutive) to consecutive internal
%   bus numbers which start at 1.
%   May be called as [i2e, bus, gen, branch] = ext2int(bus, gen, branch) if
%   area data is not available/needed.

%   MATPOWER Version 2.0
%   by Ray Zimmerman, PSERC Cornell    9/19/97
%   Copyright (c) 1997 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%% define names for columns to data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
	GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
	RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch;
[AREA_I, PRICE_REF_BUS] = idx_area;

i2e	= bus(:, BUS_I);
e2i = zeros(max(i2e), 1);
e2i(i2e) = [1:size(bus, 1)]';

bus(:, BUS_I)				= e2i( bus(:, BUS_I)			);
gen(:, GEN_BUS)				= e2i( gen(:, GEN_BUS)			);
branch(:, F_BUS)			= e2i( branch(:, F_BUS)			);
branch(:, T_BUS)			= e2i( branch(:, T_BUS)			);
if nargin > 3 & nargout > 4
	area(:, PRICE_REF_BUS)	= e2i( area(:, PRICE_REF_BUS)	);
end

return;
