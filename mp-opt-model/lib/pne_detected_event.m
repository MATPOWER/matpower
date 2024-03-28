function [ev, i] = pne_detected_event(event_log, name, zero)
% pne_detected_event - Returns the detected event of a particular name.
% ::
%
%   EV = PNE_DETECTED_EVENT(EVENT_LOG, NAME)
%   EV = PNE_DETECTED_EVENT(EVENT_LOG, NAME, ZERO)
%   [EV, I] = PNE_DETECTED_EVENT(...)
%
%   Searches through the log of detected events, returning the event with
%   the specified name, if one exists, and optionally only if it is refers
%   to a detected event zero (not an event interval).
%
%   Inputs:
%       EVENT_LOG : struct array of detected events
%       NAME : name of event of interest
%       ZERO : (optional, default = 0)
%           0 - return the event if a zero or an interval are detected for the
%               corresponding event function
%           1 - return the event only if a zero is deteced for the
%               corresponding event function
%
%   Outputs:
%       EV : entry of the EVENT_LOG corresponding to the detected event with
%            the specified NAME, or empty if NAME not in the EVENT_LOG of
%            detected events
%       I : index of EV in EVENT_LOG
%
% See pne_detect_events for details on the fields in the ``event_log`` and
% ``ev`` structs.

%   MP-Opt-Model
%   Copyright (c) 2016-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin < 3
    zero = 0;
end
ev = [];
i = 0;
for k = 1:length(event_log)
    if strcmp(event_log(k).name, name) && (event_log(k).zero || ~zero)
        ev = event_log(k);
        i = k;
    end
    break
end
