function [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch
%IDX_BRCH   Defines constants for named column indices to branch matrix.
%   [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
%   RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch
%
%   Some examples of usage, after defining the constants using the line above,
%   are:
%
%    branch(4, BR_STATUS) = 0;              % take branch 4 out of service
%    Ploss = branch(:, PF) + branch(:, PT); % compute real power loss vector
% 
%   The index, name and meaning of each column of the branch matrix is given
%   below:
%
%   columns 1-11 must be included in input matrix (in case file)
%    1  BUS_I       bus number (1 to 29997)
%    1  F_BUS       f, from bus number
%    2  T_BUS       t, to bus number
%    3  BR_R        r, resistance (p.u.)
%    4  BR_X        x, reactance (p.u.)
%    5  BR_B        b, total line charging susceptance (p.u.)
%    6  RATE_A      rateA, MVA rating A (long term rating)
%    7  RATE_B      rateB, MVA rating B (short term rating)
%    8  RATE_C      rateC, MVA rating C (emergency rating)
%    9  TAP         ratio, transformer off nominal turns ratio
%    10 SHIFT       angle, transformer phase shift angle (degrees)
%    11 BR_STATUS   initial branch status, 1 - in service, 0 - out of service
%
%   columns 12-15 are added to matrix after power flow or OPF solution
%   they are typically not present in the input matrix
%    12 PF          real power injected at "from" bus end (MW)
%    13 QF          reactive power injected at "from" bus end (MVAr)
%    14 PT          real power injected at "to" bus end (MW)
%    15 QT          reactive power injected at "to" bus end (MVAr)
%
%   columns 14-17 are added to matrix after OPF solution
%   they are typically not present in the input matrix
%                   (assume OPF objective function has units, u)
%    16 MU_SF       Kuhn-Tucker multiplier on MVA limit at "from" bus (u/MVA)
%    17 MU_ST       Kuhn-Tucker multiplier on MVA limit at "to" bus (u/MVA)

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
SHIFT       = 10;   %% angle, transformer phase shift angle (degrees)
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
