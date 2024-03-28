function efv = pne_event_target_lam(cx, opt)
% pne_event_target_lam - Event function to detect a target lambda value.
% ::
%
%   EFV = PNE_EVENT_TARGET_LAM(CX, OPT)
%
%   PNES_MASTER event function to detect the completion of the continuation
%   curve or another target value of lambda provided in OPT.stop_at.
%
%   Set OPT.stop_at a positive numerical value or to 'FULL' to use this
%   function to trigger termination of the continuation at the desired
%   value, or to trace the full curve back to lambda = 0, respectively.
%
%   Inputs:
%       CX : struct containing info about current point (continuation soln)
%       OPT - PNES_MASTER options struct
%
%   Outputs:
%       EFV : event function value
%
% See also pnes_master, pne_register_events, pne_event_nose.

%   MP-Opt-Model
%   Copyright (c) 2016-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% event function value is a scalar equal to: current lambda - target lambda
target = opt.stop_at;
if ischar(target)       %% 'FULL' trace requested
    if cx.z(end) >= 0   %% before the nose point ...
        target = -1;    %% prevent early termination (e.g. itr 1 rollback to 0)
    else                %% after the nose point ...
        target = 0;     %% terminate at lambda = 0
    end
end
efv = cx.x(end) - target;
