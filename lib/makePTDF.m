function H = makePTDF(baseMVA, bus, branch, slack)
%MAKEPTDF   Builds the DC PTDF matrix for a given choice of slack.
%   H = MAKEPTDF(MPC)
%   H = MAKEPTDF(MPC, SLACK)
%   H = MAKEPTDF(BASEMVA, BUS, BRANCH)
%   H = MAKEPTDF(BASEMVA, BUS, BRANCH, SLACK)
%
%   Returns the DC PTDF matrix for a given choice of slack. The matrix is
%   nbr x nb, where nbr is the number of branches and nb is the number of
%   buses. The SLACK can be a scalar (single slack bus) or an nb x 1 column
%   vector of weights specifying the proportion of the slack taken up at each
%   bus. If the SLACK is not specified the reference bus is used by default.
%   Bus numbers must be consecutive beginning at 1 (i.e. internal ordering).
%
%   Examples:
%       H = makePTDF(mpc);
%       H = makePTDF(baseMVA, bus, branch, 1);
%       slack = rand(size(bus, 1), 1);
%       H = makePTDF(mpc, slack);
%
%   See also MAKELODF.

%   For convenience, SLACK can also be an nb x nb matrix, where each
%   column specifies how the slack should be handled for injections
%   at that bus.

%   MATPOWER
%   Copyright (c) 2006-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% extract from MPC if necessary
if nargin < 3
    mpc     = baseMVA;
    if nargin == 2
        slack = bus;
    end
    baseMVA = mpc.baseMVA;
    bus     = mpc.bus;
    branch  = mpc.branch;
end

%% define named indices into bus matrix
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

%% use reference bus for slack by default
if nargin == 1 || nargin == 3
    slack = find(bus(:, BUS_TYPE) == REF);
    slack = slack(1);
end

%% set the slack bus to be used to compute initial PTDF
if length(slack) == 1
    slack_bus = slack;
else
    slack_bus = 1;      %% use bus 1 for temp slack bus
end

nb = size(bus, 1);
nbr = size(branch, 1);
noref   = (2:nb)';      %% use bus 1 for voltage angle reference
noslack = find((1:nb)' ~= slack_bus);

%% check that bus numbers are equal to indices to bus (one set of bus numbers)
if any(bus(:, BUS_I) ~= (1:nb)')
    error('makePTDF: buses must be numbered consecutively in bus matrix; use ext2int() to convert to internal ordering')
end

%% compute PTDF for single slack_bus
[Bbus, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);
H = zeros(nbr, nb);
H(:, noslack) = full(Bf(:, noref) / Bbus(noslack, noref));
        %%    = full(Bf(:, noref) * inv(Bbus(noslack, noref)));

%% distribute slack, if requested
if length(slack) ~= 1
    if size(slack, 2) == 1  %% slack is a vector of weights
        slack = slack/sum(slack);   %% normalize weights
        
        %% conceptually, we want to do ...
        %%    H = H * (eye(nb,nb) - slack * ones(1, nb));
        %% ... we just do it more efficiently
        v = H * slack;
        for k = 1:nb
            H(:, k) = H(:, k) - v;
        end
    else
        H = H * slack;
    end
end
