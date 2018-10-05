function delta = calc_branch_angle(mpc)
%CALC_BRANCH_ANGLE  Calculate branch angle differences across active branches
%   DELTA = CALC_BRANCH_ANGLE(MPC)
%
%   Calculates the angle difference (in degrees) across all active branches
%   in the MATPOWER case. Angles are calculated as the difference between
%   the FROM bus and the TO bus.
%
%   Input:
%       MPC   - MATPOWER case struct (can have external bus numbering)
%
%   Output:
%       DELTA - nl x 1 vector of branch angle differences Af - At, where
%               Af and At are vectors of voltage angles at "from" and "to"
%               ends of each line respectively. DELTA is 0 for out-of-service
%               branches.
%
%   See also TOGGLE_SOFTLIMS.

%   MATPOWER
%   Copyright (c) 2018, Power Systems Engineering Research Center (PSERC)
%   by Eran Schweitzer, Arizona State University
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

status = mpc.branch(:, BR_STATUS);
nl  = size(mpc.branch, 1);
nb  = size(mpc.bus, 1);
max_bus_num = max(mpc.bus(:, BUS_I));
e2i = sparse(mpc.bus(:, BUS_I), 1, 1:nb, max_bus_num, 1);   %% ext to int bus number map
bf  = full(e2i(mpc.branch(:, F_BUS)));              %% "from" bus indices
bt  = full(e2i(mpc.branch(:, T_BUS)));              %% "to" bus indices

A = sparse([1:nl,1:nl]', [bf; bt], [status; -status], nl, nb);
delta = A * mpc.bus(:, VA);                         %% angle differences
