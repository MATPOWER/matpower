function H = makePTDF(baseMVA, bus, branch, slack, bus_idx)
% makeLODF - Builds the DC PTDF matrix for a given choice of slack.
% ::
%
%   H = MAKEPTDF(MPC)
%   H = MAKEPTDF(MPC, SLACK)
%   H = MAKEPTDF(MPC, SLACK, TXFR)
%   H = MAKEPTDF(MPC, SLACK, BUS_IDX)
%   H = MAKEPTDF(BASEMVA, BUS, BRANCH)
%   H = MAKEPTDF(BASEMVA, BUS, BRANCH, SLACK)
%   H = MAKEPTDF(BASEMVA, BUS, BRANCH, SLACK, TXFR)
%   H = MAKEPTDF(BASEMVA, BUS, BRANCH, SLACK, BUS_IDX)
%
%   Returns the DC PTDF matrix for a given choice of slack. The matrix is
%   nbr x nb, where nbr is the number of branches and nb is the number of
%   buses. The SLACK can be a scalar (single slack bus) or an nb x 1 column
%   vector of weights specifying the proportion of the slack taken up at each
%   bus. If the SLACK is not specified the reference bus is used by default.
%   Bus numbers must be consecutive beginning at 1 (i.e. internal ordering).
%
%   For convenience, SLACK can also be an nb x nb matrix, where each
%   column specifies how the slack should be handled for injections
%   at that bus. This option only applies when computing the full
%   PTDF matrix (i.e. when TXFR and BUS_IDX are not provided.)
%
%   If TXFR is supplied it must be a matrix (or vector) with nb rows whose
%   columns each sum to zero, where each column defines a specific (slack
%   independent) transfer. E.g. if k-th transfer is from bus i to bus j,
%   TXFR(i, k) = 1 and TXFR(j, k) = -1. In this case H has the same number
%   of columns as TXFR.
%
%   If BUS_IDX is supplied, it contains a column vector of bus indices.
%   The columns of H correspond to these indices, but they are computed
%   individually rather than computing the full PTDF matrix and selecting
%   the desired columns.
%
%   Examples:
%       H = makePTDF(mpc);
%       H = makePTDF(baseMVA, bus, branch, 1);
%       slack = rand(size(bus, 1), 1);
%       H = makePTDF(mpc, slack);
%
%       % for transfer from bus i to bus j
%       txfr = zeros(nb, 1); txfr(i) = 1; txfr(j) = -1;
%       H = makePTDF(mpc, slack, txfr);
%
%       % for buses i and j only
%       H = makePTDF(mpc, slack, [i;j]);
%
% See also makeLODF.

%   MATPOWER
%   Copyright (c) 2006-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% extract from MPC if necessary
if isstruct(baseMVA)
    mpc = baseMVA;
    if nargin < 3
        bus_idx = [];
        if nargin < 2
            slack = [];
        else
            slack = bus;
        end
    else
        slack = bus;
        bus_idx = branch;
    end
    baseMVA = mpc.baseMVA;
    bus     = mpc.bus;
    branch  = mpc.branch;
else
    if nargin < 5
        bus_idx = [];
        if nargin < 4
            slack = [];
        end
    end
end

%% define named indices into bus matrix
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

%% use reference bus for slack by default
if isempty(slack)
    slack = find(bus(:, BUS_TYPE) == REF);
    slack = slack(1);
end

%% compute full PTDF?
nb = size(bus, 1);
nbr = size(branch, 1);
txfr = 0;               %% default assumes not just for specific transfers
if ~isempty(bus_idx)
    compute_full_H = 0;
    if size(bus_idx, 1) == nb && sum(sum(bus_idx)) < 1e-12
        txfr = 1;       %% cols of H for specific (slack independent) transfers
        dP = bus_idx;
    end
else
    compute_full_H = 1;
end

%% set the slack bus to be used to compute initial PTDF
if length(slack) == 1
    slack_bus = slack;
else
    slack_bus = 1;      %% use bus 1 for temp slack bus
end

noref   = (2:nb)';      %% use bus 1 for voltage angle reference
noslack = find((1:nb)' ~= slack_bus);

%% check that bus numbers are equal to indices to bus (one set of bus numbers)
if any(bus(:, BUS_I) ~= (1:nb)')
    error('makePTDF: buses must be numbered consecutively in bus matrix; use ext2int() to convert to internal ordering')
end

%%-----  compute PTDF for single slack_bus  -----
[Bbus, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);

%% set up injections/transfers
if compute_full_H   %% full H for all columns
%     %% old (slower) method
%     H = zeros(nbr, nb);
%     H(:, noslack) = full(Bf(:, noref) / Bbus(noslack, noref));
%             %%    = full(Bf(:, noref) * inv(Bbus(noslack, noref)));
    nbi = nb;       %% number of buses of interest (i.e. all of them)
    dP = speye(nb, nb);
else                %% compute only for specific columns
    if txfr             %% for specific transfers
        nbi = size(dP, 2);          %% number of transfers
    else                %% for subset of columns of H for specific buses
        %% if necessary add missing slacks to bus index list
        nbi0 = length(bus_idx);     %% number of original buses of interest
        bidx = bus_idx;
        if length(slack) ~= 1 && size(slack, 2) == 1
            slacks = find(slack);   %% get all slack buses
            k = find(~ismember(slacks, bus_idx));   %% find slacks not in bus_idx
            if ~isempty(k)
                bidx = [bus_idx; slacks(k)];
            end
        end
        nbi = size(bidx, 1);        %% number of buses of interest

        %% define the equivalent transfer, each column is a single transfer
        dP = accumarray([bidx (1:nbi)'], 1, [nb, nbi]);
    end
end

%% solve for change in voltage angles
dTheta = zeros(nb, nbi);
dTheta(noref, :) = Bbus(noslack, noref) \ dP(noslack, :);

%% compute corresponding change in branch flows
H = Bf * dTheta;

%%-----  distribute slack, if requested  -----
if length(slack) ~= 1 && ~txfr
    slack = slack ./ sum(slack);   %% normalize weights
    if size(slack, 2) == 1  %% slack is a vector of weights
        
        %% conceptually, we want to do ...
        %%    H = H * (eye(nb,nb) - slack * ones(1, nb));
        %% ... we just do it more efficiently
        if compute_full_H
            v = H * slack;
            for k = 1:nb
                H(:, k) = H(:, k) - v;
            end
        else
            v = H * slack(bidx);
            for k = 1:nbi
                H(:, k) = H(:, k) - v;
            end
            H = H(:, 1:nbi0);   %% remove temp cols added for missing slacks
        end
    else
        H = H - H * slack;
    end
end
