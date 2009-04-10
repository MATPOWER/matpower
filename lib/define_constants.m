%DEFINE_CONSTANTS  Defines constants for named column indices to data matrices.
%
%   DEFINE_CONSTANTS is a script file that calls IDX_AREA, IDX_BRCH,
%   IDX_BUS, IDX_COST, and IDX_GEN to define constants, especially
%   those for the named column indices for the corresponding data matrices.
%
%   Area:
%      AREA_I, PRICE_REF_BUS
%
%   Branch:
%      F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C,
%      TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST,
%      ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX
%
%   Bus:
%      PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,
%      VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN
%
%   Cost: 
%      PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST
%
%   Gen:
%      GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN,
%      MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX,
%      QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF
%
%   See also CASEFORMAT, IDX_AREA, IDX_BRCH, IDX_BUS, IDX_COST, IDX_GEN.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[AREA_I, PRICE_REF_BUS] = idx_area;

[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
