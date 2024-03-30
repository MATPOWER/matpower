function [islands, bridges, nonbridges] = find_bridges(mpc)
% find_bridges - Finds bridges in a network.
% ::
%
%   [ISLANDS, BRIDGES, NONBRIDGES] = FIND_BRIDGES(MPC)
%
%   Returns the islands, bridges and non-bridges in a network.
%   Bridges are filtered out using Tarjan's algorithm. A BRIDGE is a branch
%   whose removal breaks the island to multiple parts. The return value BRIDGES
%   is a cell array of vectors of the bus indices for each island.

%   MATPOWER
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Liangyu Zhang, Zhejiang University
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% find bridges in islands
nb  = size(mpc.bus, 1);         %% number of buses
nl  = size(mpc.branch, 1);      %% number of branches
maxb = max(mpc.bus(:, BUS_I));  %% max id of buses

% construct adjacency list of branches
e2i = sparse(mpc.bus(:, BUS_I), ones(nb, 1), 1:nb, maxb, 1);
C_on = sparse( ...
    e2i(mpc.branch(:, F_BUS)), e2i(mpc.branch(:, T_BUS)), ...
    mpc.branch(:, BR_STATUS)' .* (1:nl), ...
    nl, nl);
C_on = C_on + C_on';

% find islands
[islands, ~] = find_islands(mpc);
bridges = cell(size(islands));
nonbridges = cell(size(islands));

% find bridges per island
for i = 1:length(islands)
    island = islands{i};
    br = sparse(nl, 1, 0);

    n = 0;
    dfn = sparse(nb, 1, 0); %% order of DFS for each node
    low = sparse(nb, 1, 0); %% lowest order of node's neighbor except parent
    tarjan(island(1), -1);

    bridges{i} = find(br == 2);
    nonbridges{i} = find(br == 1);
end

%% find bridge via Tarjan's algorithm
    function tarjan(u, p)
        n = n + 1;
        dfn(u) = n;
        low(u) = n;
        v_lst = find(C_on(u, :) ~= 0);

        for ii = 1:length(v_lst)
            v = v_lst(ii);
            % ignore parent
            if (v == p)
                continue;
            end

            uv = sort([u, v]); %% sort u, v to avoid duplicate elements
            br(C_on(uv(1), uv(2))) = max(br(C_on(uv(1), uv(2))), 1); %% save edge

            if (dfn(v) == 0)
                tarjan(v, u);
                low(u) = min([low(v), low(u)]);

                if (low(v) > dfn(u))
                    % found bridge
                    br(C_on(uv(1), uv(2))) = 2;
                end
            else
                low(u) = min(dfn(v), low(u));
            end
        end
    end

end
