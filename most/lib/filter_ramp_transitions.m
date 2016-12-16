function md = filter_ramp_transitions(md0, threshold)
%FILTER_RAMP_TRANSITIONS  Creates a transition mask for ramp reserves.
%   MD = FILTER_RAMP_TRANSITIONS(MD0, THRESHOLD)
%
%   Creates a transition mask for ramping reserves, including only transitions
%   with probabilities greater than or equal to a given THRESHOLD value,
%   where the probability of the transition from j1 to j2 is taken to be the
%   conditional probability in TransMat multiplied by the conditional
%   probability of being in state j1, given that you've made it to period t.

%   MOST
%   Copyright (c) 2012-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

%% dimensions
md = md0;
nt = length(md.tstep);

%% create mask
for t = 1:nt
    if t == 1
        prob = md.tstep(t).TransMat;
    else
        prob = md.tstep(t).TransMat * prob;
    end
    md.tstep(t).TransMask = diag(prob) *  md.tstep(t).TransMat >= threshold;
end
