function [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus
%IDX_BUS   Defines variables for column indices to bus.
%   [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
%   VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus

%   MATPOWER Version 2.0
%   by Ray Zimmerman, PSERC Cornell    9/19/97
%   Copyright (c) 1996, 1997 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%% define bus types
PQ		= 1;
PV		= 2;
REF		= 3;
NONE	= 4;

%% define the indices
BUS_I		= 1;	%% bus number (1 to 29997)
BUS_TYPE	= 2;	%% bus type (1 - PQ bus, 2 - PV bus, 3 - reference bus, 4 - isolated bus)
PD			= 3;	%% Pd, real power demand (MW)
QD			= 4;	%% Qd, reactive power demand (MVAR)
GS			= 5;	%% Gs, shunt conductance (MW at V = 1.0 p.u.)
BS			= 6;	%% Bs, shunt susceptance (MVAR at V = 1.0 p.u.)
BUS_AREA	= 7;	%% area number, 1-100
VM			= 8;	%% Vm, voltage magnitude (p.u.)
VA			= 9;	%% Va, voltage angle (degrees)
BASE_KV		= 10;	%% baseKV, base voltage (kV)
ZONE		= 11;	%% zone, loss zone (1-999)
VMAX		= 12;	%% maxVm, maximum voltage magnitude (p.u.)		(not in PTI format)
VMIN		= 13;	%% minVm, minimum voltage magnitude (p.u.)		(not in PTI format)

%% included in opf solution, not necessarily in input
%% assume objective function has units, u
LAM_P		= 14;	%% Lagrange multiplier on real power mismatch (u/MW)
LAM_Q		= 15;	%% Lagrange multiplier on reactive power mismatch (u/MVAR)
MU_VMAX		= 16;	%% Kuhn-Tucker multiplier on upper voltage limit (u/p.u.)
MU_VMIN		= 17;	%% Kuhn-Tucker multiplier on lower voltage limit (u/p.u.)

return;
