function [bus, gen, branch, areas] = int2ext(i2e, bus, gen, branch, areas)
%INT2EXT   Converts internal to external bus numbering.
%   [bus, gen, branch, areas] = int2ext(i2e, bus, gen, branch, areas) converts
%   from the consecutive internal bus numbers back to the originals.
%   May be called as [bus, gen, branch] = int2ext(i2e, bus, gen, branch) if
%   area data is not available/needed.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define names for columns to data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[AREA_I, PRICE_REF_BUS] = idx_area;

bus(:, BUS_I)               = i2e( bus(:, BUS_I)            );
gen(:, GEN_BUS)             = i2e( gen(:, GEN_BUS)          );
branch(:, F_BUS)            = i2e( branch(:, F_BUS)         );
branch(:, T_BUS)            = i2e( branch(:, T_BUS)         );
if nargin > 4 & nargout > 3
    areas(:, PRICE_REF_BUS) = i2e( areas(:, PRICE_REF_BUS)  );
end

return;
