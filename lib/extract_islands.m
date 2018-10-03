function mpck = extract_islands(mpc, varargin)
%function mpck = extract_islands(mpc, groups, k, custom)
%EXTRACT_ISLANDS Extracts each island in a network with islands
%   MPC_ARRAY = EXTRACT_ISLANDS(MPC)
%   MPC_ARRAY = EXTRACT_ISLANDS(MPC, GROUPS)
%   MPC_K = EXTRACT_ISLANDS(MPC, K)
%   MPC_K = EXTRACT_ISLANDS(MPC, GROUPS, K)
%   MPC_K = EXTRACT_ISLANDS(MPC, K, CUSTOM)
%   MPC_K = EXTRACT_ISLANDS(MPC, GROUPS, K, CUSTOM)
%
%   Returns a cell array of MATPOWER case structs for each island in
%   the input case struct. If the optional second argument is a cell
%   array GROUPS it is assumed to be a cell array of vectors of bus
%   indices for each island (as returned by FIND_ISLANDS). Providing
%   the GROUPS avoids the need for another traversal of the network
%   connectivity and can save a significant amount of time on very
%   large systems. If an additional argument K is included, it indicates
%   which island(s) to return and the return value is a single case
%   struct, rather than a cell array. If K is a scalar or vector, it
%   it specifies the index(indices) of the island(s) to include in
%   the resulting case file. K can also be the string 'all' which
%   will include all islands. This is the same as simply eliminating
%   all isolated buses.
%
%   A final optional argument CUSTOM is a struct that can be used to
%   indicate custom fields of MPC from which to extract data
%   corresponding to buses generators, branches or DC lines. It has
%   the following structure:
%
%       CUSTOM.<ORDERING>{DIM} = FIELDS
%
%   <ORDERING> is either 'bus', 'gen', 'branch' or 'dcline' and
%   indicates that dimension DIM of FIELDS has dimensions
%   corresponding to this <ORDERING> and should have the appropriate
%   dimension extracted as well. FIELDS is a cell array, where
%   each element is either a single string (field name of MPC) or
%   a cell array of strings (nested fields of MPC).
%
%   Examples:
%       Extract each island into it's own case struct:
%           mpc_list = extract_islands(mpc);
%
%       Extract the 2nd (that is, 2nd largest) island:
%           mpc2 = extract_islands(mpc, 2);
%
%       Extract the first and 3rd islands without a re-traverals of the
%       network:
%           groups = find_islands(mpc);
%           mpc1 = extract_islands(mpc, groups, 1);
%           mpc3 = extract_islands(mpc, groups, 3);
%
%       Extract the 2nd island, including custom fields, where
%       mpc.bus_label{b} contains a label for bus b, and mpc.gen_name{g},
%       mpc.emissions.rate(g, :), and mpc.genloc(:, g) contain,
%       respectively, the generator's name, emission rates and
%       location coordinates:
%           custom.bus{1} = {'bus_label'};
%           custom.gen{1} = {'gen_name', {'emissions', 'rate'}};
%           custom.gen{2} = {'genloc'};
%           mpc = extract_islands(mpc, 1, custom);
%
%       Note: Fields bus_name, gentype and genfuel are handled automatically
%             and do not need to be included in custom.
%
%   See also FIND_ISLANDS, CASE_INFO, CONNECTED_COMPONENTS.

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
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
c = idx_dcline;

%% set up connectivity matrices
nb  = size(mpc.bus, 1);     %% number of buses
nl  = size(mpc.branch, 1);  %% number of branches
if isfield(mpc, 'dcline')   %% number of DC lines
    ndc = size(mpc.dcline, 1);
else
    ndc = 0;
end
ng  = size(mpc.gen, 1);     %% number of dispatchable injections

e2i = sparse(mpc.bus(:, BUS_I), ones(nb, 1), 1:nb, max(mpc.bus(:, BUS_I)), 1);
C_on = sparse(1:nl, e2i(mpc.branch(:, F_BUS)), -mpc.branch(:, BR_STATUS), nl, nb) + ...
    sparse(1:nl, e2i(mpc.branch(:, T_BUS)),  mpc.branch(:, BR_STATUS), nl, nb);
C = sparse(1:nl, e2i(mpc.branch(:, F_BUS)), -1, nl, nb) + ...
    sparse(1:nl, e2i(mpc.branch(:, T_BUS)),  1, nl, nb);
if ndc
    Cdc_on = sparse(1:ndc, e2i(mpc.dcline(:, c.F_BUS)), -mpc.dcline(:, c.BR_STATUS), ndc, nb) + ...
        sparse(1:ndc, e2i(mpc.dcline(:, c.T_BUS)),  mpc.dcline(:, c.BR_STATUS), ndc, nb);
    Cdc = sparse(1:ndc, e2i(mpc.dcline(:, c.F_BUS)), -1, ndc, nb) + ...
        sparse(1:ndc, e2i(mpc.dcline(:, c.T_BUS)),  1, ndc, nb);
end
Cg_on = sparse(1:ng, e2i(mpc.gen(:, GEN_BUS)), mpc.gen(:, GEN_STATUS), ng, nb);
Cg = sparse(1:ng, e2i(mpc.gen(:, GEN_BUS)), 1, ng, nb);

if nnz(C)
    n = length(varargin);
    if n >= 1 && iscell(varargin{1})
        groups = varargin{1};
        z = 1;
    else
        groups = {};
        z = 0;
    end
    if z+1 <= n
        k = varargin{z+1};
    else
        k = [];
    end
    if z+2 <= n
        custom = varargin{z+2};
    else
        custom = struct();
    end
    
    %% find islands, if not provided
    if isempty(groups)
        groups = connected_components(C_on);
    end
    
    %% check inputs
    if isempty(k)
        g1 = 1;
        gn = length(groups);
    else
        if ischar(k)
            if strcmp(upper(k), 'ALL')
                k = (1:length(groups))';
            else
                error('extract_islands: K = ''%s'' is not a valid input', k);
            end
        end
        if max(k) > length(groups)
            error('extract_islands: cannot extract island %d, network has only %d islands', ...
                    max(k), length(groups));
        end
        if length(k) > 1        %% extract multiple islands as one case
            tmpgroup = groups{k(1)};
            for j = 2:length(k)
                tmpgroup = union(tmpgroup, groups{k(j)});
            end
            groups = { tmpgroup };
            g1 = 1;
            gn = 1;
        else                    %% extract single island
            g1 = k;
            gn = k;
        end
    end
    
    %% extract islands
    for i = g1:gn
        kk  = i-g1+1;
        b  = groups{i};                     %% buses in group i
        %% branches with both ends in group i
        ibr = find(sum(abs(C(:, b)), 2) & ~sum(C(:, b), 2));
        ig  = find(sum(Cg(:, b), 2));       %% gens in group i
        %% DC lines with both ends in group i
        if ndc
            idc = find(sum(abs(Cdc(:, b)), 2) & ~sum(Cdc(:, b), 2));
        else
            idc = [];
        end
    
        mpck{kk}        = mpc;
        mpck{kk}.bus    = mpc.bus(b, :);
        mpck{kk}.branch = mpc.branch(ibr, :);
        mpck{kk}.gen    = mpc.gen(ig, :);
        if isfield(mpck{kk}, 'gencost')
            if size(mpck{kk}.gencost, 1) == 2*ng
                mpck{kk}.gencost = mpc.gencost([ig; ng+ig], :);
            else
                mpck{kk}.gencost = mpc.gencost(ig, :);
            end
        end
        if isfield(mpck{kk}, 'gentype')
            mpck{kk}.gentype = mpc.gentype(ig);
        end
        if isfield(mpck{kk}, 'genfuel')
            mpck{kk}.genfuel = mpc.genfuel(ig);
        end
        if isfield(mpck{kk}, 'bus_name')
            mpck{kk}.bus_name = mpc.bus_name(b);
        end
        if ndc
            mpck{kk}.dcline = mpc.dcline(idc, :);
            if isfield(mpck{kk}, 'dclinecost')
                mpck{kk}.dclinecost = mpc.dclinecost(idc, :);
            end
        end

        %% handle custom fields
        orderings = {'bus', 'gen', 'branch', 'dcline'};
        indices = {b, ig, ibr, idc};

        for n = 1:length(orderings)
            ord = orderings{n};
            if isfield(custom, ord)
                for dim = 1:length(custom.(ord))
                    for j = 1:length(custom.(ord){dim})
                        s = [];
                        field = custom.(ord){dim}{j};
                        if ischar(field)
                            field = { field };
                        end

                        tmp = mpck{kk}; %% check this for presence of sub-fields
                        skip = 0;
                        for i = 1:length(field)
                            s(i).type = '.';
                            s(i).subs = field{i};
                            if isfield(tmp, field{i}) && ~isempty(tmp.(field{i}))
                                %% have sub-field, continue
                                tmp = tmp.(field{i});
                            else
                                %% sub-field doesn't exist, skip it
                                skip = 1;
                                break;
                            end
                        end
                        if ~skip
                            mpck{kk} = subsasgn(mpck{kk}, s, get_reorder(subsref(mpck{kk}, s), indices{n}, dim));
                        end
                    end
                end
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
