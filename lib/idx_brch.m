function [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch
%IDX_BRCH   Defines variables for column indices to branch.
%   [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
%   RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define the indices
F_BUS       = 1;    %% f, from bus number
T_BUS       = 2;    %% t, to bus number
BR_R        = 3;    %% r, resistance (p.u.)
BR_X        = 4;    %% x, reactance (p.u.)
BR_B        = 5;    %% b, total line charging susceptance (p.u.)
RATE_A      = 6;    %% rateA, MVA rating A (long term rating)
RATE_B      = 7;    %% rateB, MVA rating B (short term rating)
RATE_C      = 8;    %% rateC, MVA rating C (emergency rating)
TAP         = 9;    %% ratio, transformer off nominal turns ratio
SHIFT       = 10;   %% angle, transformer phase shift angle
BR_STATUS   = 11;   %% initial branch status, 1 - in service, 0 - out of service

%% included in power flow solution, not necessarily in input
PF          = 12;   %% real power injected at "from" bus end (MW)       (not in PTI format)
QF          = 13;   %% reactive power injected at "from" bus end (MVAr) (not in PTI format)
PT          = 14;   %% real power injected at "to" bus end (MW)         (not in PTI format)
QT          = 15;   %% reactive power injected at "to" bus end (MVAr)   (not in PTI format)

%% included in opf solution, not necessarily in input
%% assume objective function has units, u
MU_SF       = 16;   %% Kuhn-Tucker multiplier on MVA limit at "from" bus (u/MVA)
MU_ST       = 17;   %% Kuhn-Tucker multiplier on MVA limit at "to" bus (u/MVA)

return;
