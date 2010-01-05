%DEFINE_CONSTANTS  Defines constants for named column indices to data matrices.
%
%   define_constants
%
%   This is simply a convenience script that defines the constants
%   listed below, consisting primarily of named indices for the
%   columns of the data matrices: bus, branch, gen and gencost.
%   This includes input columns defined in caseformat as well as
%   columns that are added in the power flow and OPF output.
%
%   bus:
%      PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,
%      VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN
%
%   branch:
%      F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C,
%      TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST,
%      ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX
%
%   gen:
%      GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN,
%      MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX,
%      QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF
%
%   gencost: 
%      PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST
%
%   See the help for caseformat, idx_bus, idx_brch, idx_gen and
%   idx_cost for details on the meaning of these constants. Internally
%   define_constants calls idx_bus, idx_brch, idx_gen and idx_cost.
%   In performance sensitive code, such as internal MATPOWER functions
%   that are called frequently, it is preferred to call these
%   functions directly rather than using the define_constants script,
%   which is less efficient.
%
%   This script is included for convenience for interactive use or
%   for high-level code where maximum performance is not a concern.

%   MATPOWER
%   $Id$
%   by Doug Mitarotonda & Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009-2010 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
%[AREA_I, PRICE_REF_BUS] = idx_area;    %% deprecated
