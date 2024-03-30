function LODF = makeLODF(branch, PTDF, mask_bridge)
% makeLODF - Builds the line outage distribution factor matrix.
% ::
%
%   LODF = MAKELODF(BRANCH, PTDF)
%   LODF = MAKELODF(MPC, PTDF)
%   LODF = MAKELODF(MPC, PTDF, MASK_BRIDGE)
% 
%   Returns the DC model line outage distribution factor matrix corresponding
%   to a given PTDF matrix. The LODF matrix is nbr x nbr, where nbr is the
%   number of branches. If the optional MASK_BRIDGE argument is true, columns
%   corresponding to bridge branches (those whose removal result in
%   islanding) are replaced with NaN.
%
%   Example:
%       H = makePTDF(mpc);
%       LODF = makeLODF(branch, H);
%       LODF = makeLODF(mpc, H);
%
%       % mask bridge branches in LODF
%       makeLODF(mpc, H, 1);
%
% See also makePTDF, find_bridges.

%   MATPOWER
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell,
%   and Liangyu Zhang, Zhejiang University
%
%   Modified by Liangyu Zhang
%   2023.05.19 Mask corresponding columns of network bridge branches

%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if isstruct(branch)
    mpc = branch;
    branch = mpc.branch;
end

if nargin < 3
    mask_bridge = 0;
elseif mask_bridge && ~exist('mpc', 'var')
    error('MPC is required when mask_bridge is enabled.');
end

[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

[nl, nb] = size(PTDF);
f = branch(:, F_BUS);
t = branch(:, T_BUS);
Cft =  sparse([f; t], [1:nl 1:nl]', [ones(nl, 1); -ones(nl, 1)], nb, nl);

H = PTDF * Cft;
h = diag(H, 0);
LODF = H ./ (ones(nl, nl) - ones(nl, 1) * h');
LODF = LODF - diag(diag(LODF)) - eye(nl, nl);

if mask_bridge
    [~, bridges, nonbridges] = find_bridges(mpc);
    for i = 1:length(bridges)
        br = bridges{i};
        nonbr = nonbridges{i};
        LODF(:, br) = NaN;  % mask bridge branch columns
        LODF(setdiff(1:nb, [br; nonbr]), nonbr) = NaN; % one island should not affect other islands
    end
end
