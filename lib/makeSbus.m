function Sbus = makeSbus(baseMVA, bus, gen)
%MAKESBUS   Builds the vector of complex bus power injections.
%   Sbus = makeSbus(baseMVA, bus, gen) returns the vector of complex bus
%   power injections, that is, generation minus load. Power is expressed
%   in per unit.

%   MATPOWER Version 2.5b3
%   by Ray Zimmerman, PSERC Cornell    7/15/99
%   Copyright (c) 1996-1999 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%% constants
j = sqrt(-1);

%% define named indices into bus, gen matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
	GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;

%% generator info
on = find(gen(:, GEN_STATUS) > 0);		%% which generators are on?
gbus = gen(on, GEN_BUS);				%% what buses are they at?

%% form net complex bus power injection vector
nb = size(bus, 1);
ngon = size(on, 1);
Cg = sparse(gbus, [1:ngon]', ones(ngon, 1), nb, ngon);	%% connection matrix
														%% element i, j is 1 if
														%% gen on(j) at bus i is ON
Sbus =	( Cg * (gen(on, PG) + j * gen(on, QG)) ...	%% power injected by generators
			- (bus(:, PD) + j * bus(:, QD)) ) / ...	%% plus power injected by loads
		baseMVA;									%% converted to p.u.

return;
