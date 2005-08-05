function [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, QMAX2, QMIN2, ...
    RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen
%IDX_GEN   Defines constants for named column indices to gen matrix.
%
%   [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, ...
%   PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, QMAX2, QMIN2, ...
%   RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen
%
%   Some examples of usage, after defining the constants using the line above,
%   are:
%
%    Pg = gen(4, PG);   % get the real power output of generator 4
%    gen(:, PMIN) = 0;  % set to zero the minimum real power limit of all gens
% 
%   The index, name and meaning of each column of the gen matrix is given
%   below:
%
%   columns 1-03 must be included in input matrix (in case file)
%    1  GEN_BUS     bus number
%    2  PG          Pg, real power output (MW)
%    3  QG          Qg, reactive power output (MVAr)
%    4  QMAX        Qmax, maximum reactive power output (MVAr)
%    5  QMIN        Qmin, minimum reactive power output (MVAr)
%    6  VG          Vg, voltage magnitude setpoint (p.u.)
%    7  MBASE       mBase, total MVA base of machine, defaults to baseMVA
%    8  GEN_STATUS  status, 1 - in service, 0 - out of service
%    9  PMAX        Pmax, maximum real power output (MW)
%    10 PMIN        Pmin, minimum real power output (MW)
%    11 QMAX2       Qmax2, maximum reactive power output at Pmax (MVAr)
%    12 QMIN2       Qmin2, minimum reactive power output at Pmax (MVAr)
%    13 RAMP_AGC    ramp rate for load following/AGC (MW/min)
%    14 RAMP_10     ramp rate for 10 minute reserves (MW)
%    15 RAMP_30     ramp rate for 30 minute reserves (MW)
%    16 RAMP_Q      ramp rate for reactive power (2 sec timescale) (MVAr/min)
%    17 APF         area participation factor
%   
%   columns 17-20 are added to matrix after OPF solution
%   they are typically not present in the input matrix
%                   (assume OPF objective function has units, u)
%    18 MU_PMAX     Kuhn-Tucker multiplier on upper Pg limit (u/MW)
%    19 MU_PMIN     Kuhn-Tucker multiplier on lower Pg limit (u/MW)
%    20 MU_QMAX     Kuhn-Tucker multiplier on upper Qg limit (u/MVAr)
%    21 MU_QMIN     Kuhn-Tucker multiplier on lower Qg limit (u/MVAr)

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2005 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define the indices
GEN_BUS     = 1;    %% bus number
PG          = 2;    %% Pg, real power output (MW)
QG          = 3;    %% Qg, reactive power output (MVAr)
QMAX        = 4;    %% Qmax, maximum reactive power output at Pmin (MVAr)
QMIN        = 5;    %% Qmin, minimum reactive power output at Pmin (MVAr)
VG          = 6;    %% Vg, voltage magnitude setpoint (p.u.)
MBASE       = 7;    %% mBase, total MVA base of this machine, defaults to baseMVA
GEN_STATUS  = 8;    %% status, 1 - machine in service, 0 - machine out of service
PMAX        = 9;    %% Pmax, maximum real power output (MW)
PMIN        = 10;   %% Pmin, minimum real power output (MW)
QMAX2       = 11;   %% maximum reactive power output at Pmax (MVAr)
QMIN2       = 12;   %% minimum reactive power output at Pmax (MVAr)
RAMP_AGC    = 13;   %% ramp rate for load following/AGC (MW/min)
RAMP_10     = 14;   %% ramp rate for 10 minute reserves (MW)
RAMP_30     = 15;   %% ramp rate for 30 minute reserves (MW)
RAMP_Q      = 16;   %% ramp rate for reactive power (2 sec timescale) (MVAr/min)
APF         = 17;   %% area participation factor

%% included in opf solution, not necessarily in input
%% assume objective function has units, u
MU_PMAX     = 18;   %% Kuhn-Tucker multiplier on upper Pg limit (u/MW)
MU_PMIN     = 19;   %% Kuhn-Tucker multiplier on lower Pg limit (u/MW)
MU_QMAX     = 20;   %% Kuhn-Tucker multiplier on upper Qg limit (u/MVAr)
MU_QMIN     = 21;   %% Kuhn-Tucker multiplier on lower Qg limit (u/MVAr)

%% Note: When QMAX2 and QMAX are not equal (generator capability curve is
%% not simply a box) and the upper Qg limit is binding, the multiplier on
%% this constraint is split into it's P and Q components and combined with the
%% appropriate MU_Pxxx and MU_Qxxx values. Likewise for the lower Q limits.

return;
