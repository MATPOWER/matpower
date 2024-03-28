function efv = pne_event_nose(cx, opt)
% pne_event_nose - Event function to detect the limit or nose point.
% ::
%
%   EFV = PNE_EVENT_NOSE(CX, OPT)
%
%   PNE_MASTER event function to detect the limit or nose point of the
%   continuation curve, based on the sign of the lambda component of the
%   tangent vector.
%
%   Set OPT.stop_at to 'NOSE' to use this function to trigger termination
%   of the continuation at the limit or nose point.
%
%   Inputs:
%       CX : struct containing info about current point (continuation soln)
%       OPT - PNES_MASTER options struct
%
%   Outputs:
%       EFV : event function value
%
% See also pnes_master, pne_register_events, pne_event_target_lam.

%   MP-Opt-Model
%   Copyright (c) 2016-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% event function value is dlam, the last element of the
%% normalized tangent vector at the current soln
efv = cx.z(end);
