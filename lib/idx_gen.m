function [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
	MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen
%IDX_GEN   Defines variables for column indices to gen.
%   [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, ...
%   PMAX, PMIN, START, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2003 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%% define the indices
GEN_BUS		= 1;	%% bus number
PG			= 2;	%% Pg, real power output (MW)
QG			= 3;	%% Qg, reactive power output (MVAr)
QMAX		= 4;	%% Qmax, maximum reactive power output (MVAr)
QMIN		= 5;	%% Qmin, minimum reactive power output (MVAr)
VG			= 6;	%% Vg, voltage magnitude setpoint (p.u.)
MBASE		= 7;	%% mBase, total MVA base of this machine, defaults to baseMVA
GEN_STATUS	= 8;	%% status, 1 - machine in service, 0 - machine out of service
PMAX		= 9;	%% Pmax, maximum real power output (MW)
PMIN		= 10;	%% Pmin, minimum real power output (MW)

%% included in opf solution, not necessarily in input
%% assume objective function has units, u
MU_PMAX		= 11;	%% Kuhn-Tucker multiplier on upper Pg limit (u/MW)
MU_PMIN		= 12;	%% Kuhn-Tucker multiplier on lower Pg limit (u/MW)
MU_QMAX		= 13;	%% Kuhn-Tucker multiplier on upper Qg limit (u/MVAr)
MU_QMIN		= 14;	%% Kuhn-Tucker multiplier on lower Qg limit (u/MVAr)

return;
