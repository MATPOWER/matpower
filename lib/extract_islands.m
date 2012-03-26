function mpck = extract_islands(mpc, groups, k)
%EXTRACT_ISLANDS Extracts each island in a network with islands
%   MPC_ARRAY = EXTRACT_ISLANDS(MPC)
%   MPC_ARRAY = EXTRACT_ISLANDS(MPC, GROUPS)
%   MPC_K = EXTRACT_ISLANDS(MPC, K)
%   MPC_K = EXTRACT_ISLANDS(MPC, GROUPS, K)
%
%   Returns a cell array of MATPOWER case structs for each island in
%   the input case struct. If the optional second argument is a cell
%   array GROUPS it is assumed to be a cell array of vectors of bus
%   indices for each island (as returned by FIND_ISLANDS). Providing
%   the GROUPS avoids the need for another traversal of the network
%   connectivity and can save a significant amount of time on very
%   large systems. If an additional scalar argument K is included, it
%   indicates which island to return and the return value is the single
%   corresponding case struct itself.
%
%   See also FIND_ISLANDS, CONNECTED_COMPONENTS.

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
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% set up connectivity matrices
nb  = size(mpc.bus, 1);     %% number of buses
nl  = size(mpc.branch, 1);  %% number of branches
ng  = size(mpc.gen, 1);     %% number of dispatchable injections

e2i = sparse(mpc.bus(:, BUS_I), ones(nb, 1), 1:nb, max(mpc.bus(:, BUS_I)), 1);
C = sparse(1:nl, e2i(mpc.branch(:, F_BUS)), -mpc.branch(:, BR_STATUS), nl, nb) + ...
    sparse(1:nl, e2i(mpc.branch(:, T_BUS)),  mpc.branch(:, BR_STATUS), nl, nb);
Cg = sparse(1:ng, e2i(mpc.gen(:, GEN_BUS)), mpc.gen(:, GEN_STATUS), ng, nb);

if nnz(C)
    if nargin == 2
        if iscell(groups)
            k = [];
        else
            k = groups;
            groups = {};
        end
    elseif nargin == 1
        groups = {};
        k = [];
    end
    
    %% find islands, if not provided
    if isempty(groups)
        groups = connected_components(C);
    end
    
    %% check inputs
    if isempty(k)
        g1 = 1;
        gn = length(groups);
    else
        if k > length(groups)
            error('extract_islands: cannot extract island %d, network has only %d islands', ...
                    k, length(groups));
        end
        g1 = k;
        gn = k;
    end
    
    %% extract islands
    for i = g1:gn
        kk  = i-g1+1;
        ib  = groups{i};                    %% buses in group i
        ibr = find(sum(abs(C(:, ib)), 2));  %% branches in group i
        ig  = find(sum(Cg(:, ib), 2));      %% gens in group i
    
        mpck{kk}        = mpc;
        mpck{kk}.bus    = mpc.bus(ib, :);
        mpck{kk}.branch = mpc.branch(ibr, :);
        mpck{kk}.gen    = mpc.gen(ig, :);
        if isfield(mpck{kk}, 'gencost')
            if size(mpck{kk}.gencost, 1) == 2*ng
                mpck{kk}.gencost = mpc.gencost([ig; ng+ig], :);
            else
                mpck{kk}.gencost = mpc.gencost(ig, :);
            end
        end
    end
    
    %% convert from cell array to single MATPOWER case struct as appropriate
    if ~isempty(k)
        mpck = mpck{1};
    end
else
    mpck = [];
end