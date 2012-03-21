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
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2012 by Power System Engineering Research Center (PSERC)
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
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

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
C = sparse(1:nl, e2i(mpc.branch(:, F_BUS)), -mpc.branch(:, BR_STATUS), nl, nb) + ...
    sparse(1:nl, e2i(mpc.branch(:, T_BUS)),  mpc.branch(:, BR_STATUS), nl, nb);

if isempty(C)
    groups = [];
    isolated = 1:nb;
else
    [groups, isolated] = connected_components(C);
end
