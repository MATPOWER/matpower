function Sbus = makeSbus(baseMVA, bus, gen)
%MAKESBUS   Builds the vector of complex bus power injections.
%   Sbus = makeSbus(baseMVA, bus, gen) returns the vector of complex bus
%   power injections, that is, generation minus load. Power is expressed
%   in per unit.

%   MATPOWER Version 2.0
%   by Ray Zimmerman, PSERC Cornell    9/19/97
%   Copyright (c) 1996, 1997 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%% constants
j = sqrt(-1);

%% define named indices into bus, gen matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
	GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;

%% generator info
on = find(gen(:, GEN_STATUS));				%% which generators are on?
gbus = gen(on, GEN_BUS);					%% what buses are they at?

%% form net complex bus power injection vector
Sbus = -(bus(:, PD) + j * bus(:, QD));						%% power injected by loads
Sbus(gbus) = Sbus(gbus) + gen(on, PG) + j * gen(on, QG);	%% plus generation
Sbus = Sbus / baseMVA;										%% convert to p.u.

return;
