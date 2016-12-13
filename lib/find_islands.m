function [groups, isolated] = find_islands(mpc)
%FIND_ISLANDS  Finds islands in a network
%   GROUPS = FIND_ISLANDS(MPC)
%   [GROUPS, ISOLATED] = FIND_ISLANDS(MPC)
%
%   Returns the islands in a network. The return value GROUPS
%   is a cell array of vectors of the bus indices for each island.
%   The second and optional return value ISOLATED is a vector of
%   indices of isolated buses that have no connecting branches.
%
%   See also EXTRACT_ISLANDS, CONNECTED_COMPONENTS.

%   TODO: add handling of DC lines

%   MATPOWER
%   Copyright (c) 2012-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
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

%% find islands
nb  = size(mpc.bus, 1);     %% number of buses
nl  = size(mpc.branch, 1);  %% number of branches

e2i = sparse(mpc.bus(:, BUS_I), ones(nb, 1), 1:nb, max(mpc.bus(:, BUS_I)), 1);
C_on = sparse(1:nl, e2i(mpc.branch(:, F_BUS)), -mpc.branch(:, BR_STATUS), nl, nb) + ...
       sparse(1:nl, e2i(mpc.branch(:, T_BUS)),  mpc.branch(:, BR_STATUS), nl, nb);

if nnz(C_on)
    [groups, isolated] = connected_components(C_on);
else
    groups = [];
    isolated = 1:nb;
end
