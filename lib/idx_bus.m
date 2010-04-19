function [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus
%IDX_BUS   Defines constants for named column indices to bus matrix.
%   Example:
%
%   [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
%   VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
%
%   Some examples of usage, after defining the constants using the line above,
%   are:
%
%    Pd = bus(4, PD);       % get the real power demand at bus 4
%    bus(:, VMIN) = 0.95;   % set the min voltage magnitude to 0.95 at all buses
% 
%   The index, name and meaning of each column of the bus matrix is given
%   below:
%
%   columns 1-13 must be included in input matrix (in case file)
%    1  BUS_I       bus number (positive integer)
%    2  BUS_TYPE    bus type (1 = PQ, 2 = PV, 3 = ref, 4 = isolated)
%    3  PD          Pd, real power demand (MW)
%    4  QD          Qd, reactive power demand (MVAr)
%    5  GS          Gs, shunt conductance (MW demanded at V = 1.0 p.u.)
%    6  BS          Bs, shunt susceptance (MVAr injected at V = 1.0 p.u.)
%    7  BUS_AREA    area number, (positive integer)
%    8  VM          Vm, voltage magnitude (p.u.)
%    9  VA          Va, voltage angle (degrees)
%    10 BASE_KV     baseKV, base voltage (kV)
%    11 ZONE        zone, loss zone (positive integer)
%    12 VMAX        maxVm, maximum voltage magnitude (p.u.)
%    13 VMIN        minVm, minimum voltage magnitude (p.u.)
%   
%   columns 14-17 are added to matrix after OPF solution
%   they are typically not present in the input matrix
%                   (assume OPF objective function has units, u)
%    14 LAM_P       Lagrange multiplier on real power mismatch (u/MW)
%    15 LAM_Q       Lagrange multiplier on reactive power mismatch (u/MVAr)
%    16 MU_VMAX     Kuhn-Tucker multiplier on upper voltage limit (u/p.u.)
%    17 MU_VMIN     Kuhn-Tucker multiplier on lower voltage limit (u/p.u.)
% 
%   additional constants, used to assign/compare values in the BUS_TYPE column
%    1  PQ    PQ bus
%    2  PV    PV bus
%    3  REF   reference bus
%    4  NONE  isolated bus
%
%   See also DEFINE_CONSTANTS.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2010 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as M-files and MEX-files) available in a
%   Matlab (or compatible) environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

%% define bus types
PQ      = 1;
PV      = 2;
REF     = 3;
NONE    = 4;

%% define the indices
BUS_I       = 1;    %% bus number (1 to 29997)
BUS_TYPE    = 2;    %% bus type (1 - PQ bus, 2 - PV bus, 3 - reference bus, 4 - isolated bus)
PD          = 3;    %% Pd, real power demand (MW)
QD          = 4;    %% Qd, reactive power demand (MVAr)
GS          = 5;    %% Gs, shunt conductance (MW at V = 1.0 p.u.)
BS          = 6;    %% Bs, shunt susceptance (MVAr at V = 1.0 p.u.)
BUS_AREA    = 7;    %% area number, 1-100
VM          = 8;    %% Vm, voltage magnitude (p.u.)
VA          = 9;    %% Va, voltage angle (degrees)
BASE_KV     = 10;   %% baseKV, base voltage (kV)
ZONE        = 11;   %% zone, loss zone (1-999)
VMAX        = 12;   %% maxVm, maximum voltage magnitude (p.u.)      (not in PTI format)
VMIN        = 13;   %% minVm, minimum voltage magnitude (p.u.)      (not in PTI format)

%% included in opf solution, not necessarily in input
%% assume objective function has units, u
LAM_P       = 14;   %% Lagrange multiplier on real power mismatch (u/MW)
LAM_Q       = 15;   %% Lagrange multiplier on reactive power mismatch (u/MVAr)
MU_VMAX     = 16;   %% Kuhn-Tucker multiplier on upper voltage limit (u/p.u.)
MU_VMIN     = 17;   %% Kuhn-Tucker multiplier on lower voltage limit (u/p.u.)
