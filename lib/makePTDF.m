function H = makePTDF(baseMVA, bus, branch, slack)
%MAKEPTDFDC   Builds the DC PTDF matrix for a given choice of slack.
%   H = makePTDF(baseMVA, bus, branch, slack) returns the DC PTDF
%   matrix for a given choice of slack. The matrix is nbr x nb, where
%   nbr is the number of branches and nb is the number of buses. The slack
%   can be a scalar (single slack bus) or an nb x 1 column vector of
%   weights specifying the proportion of the slack taken up at each bus.

%   For convenience, slack can also be an nb x nb matrix, where each
%   column specifies how the slack should be handled for injections
%   at that bus.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2006 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

%% set the slack bus to be used to compute initial PTDF
if length(slack) == 1
    slack_bus = slack;
else
    slack_bus = 1;      %% use bus 1 for temp slack bus
end

nb = size(bus, 1);
nbr = size(branch, 1);
noref   = [2:nb]';      %% use bus 1 for voltage angle reference
noslack = find([1:nb]' ~= slack_bus);

%% check that bus numbers are equal to indices to bus (one set of bus numbers)
if any(bus(:, BUS_I) ~= [1:nb]')
    error('makePTDF: buses must appear in order by bus number')
end

%% compute PTDF for single slack_bus
[Bbus, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);
H = zeros(nbr, nb);
H(:, noslack) = full(Bf(:, noref) * inv(Bbus(noslack, noref)));

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

return;
