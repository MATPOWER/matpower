function reg_ev = pne_register_events(my_events, opt, reg_ev)
% pne_register_events - Register PNE event functions.
% ::
%
%   REG_EV = PNE_REGISTER_EVENTS(MY_EVENTS, OPT)
%   REG_EV = PNE_REGISTER_EVENTS(MY_EVENTS, OPT, REG_EV)
%
%   Registers event functions to be called by PNES_MASTER.
%
%   Inputs:
%       MY_EVENTS : a cell array of the form {NAME, FCN, TOL, LOCATE}, or
%           a cell array of such cell arrays, TOL and LOCATE are optional:
%           NAME - char array with a unique event name
%           FCN - function handle to the event detection function
%           TOL (OPT.default_event_tol) - optional, scalar or vector of same
%               dimension as event function return value of tolerance for
%               detecting the event, i.e. abs(val) <= tol
%           LOCATE (1) - optional, flag indicating whether the event requests
%               a rollback step to locate the event function zero
%       OPT - PNES_MASTER options struct
%       REG_EV : (optional) struct array containing existing registered event
%           functions
%
%   Outputs:
%       REG_EV : updated struct containing registered event functions, with
%           fields name, fcn, tol, locate, corresponding to respective inputs
%
%   User Defined PNES_MASTER Event Functions:
%       The user can define their own event detection functions which take
%       the same form as PNE_EVENT_TARGET_LAM. These are specified via the
%       'events' option (OPT.events) to PNES_MASTER, which takes the same
%       form as MY_EVENTS above.
%
% See also pnes_master, pne_event_nose, pne_event_target_lam.

%   MP-Opt-Model
%   Copyright (c) 2016-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% initialize registered event function struct
if nargin < 3
    reg_ev = [];
end

%% convert input to cell array of cell arrays, if necessary
if ~iscell(my_events{1})
    my_events = {my_events};
end

for k = 1:length(my_events)
    params = my_events{k};      %% event function parameters
    name = params{1};
    if length(params) > 2
        tol = params{3};
        if length(params) > 3
            locate = params{4};
        else
            locate = 1;
        end
    else
        tol = opt.default_event_tol;
    end
    e = struct( ...
            'name', name, ...
            'fcn', params{2}, ...
            'tol', tol, ...
            'locate', locate ...
        );
    if isempty(reg_ev)
        reg_ev = e;
    else
        nef = length(reg_ev);
        for k = 1:nef
            if strcmp(reg_ev(k).name, name)
                error('pne_register_events: duplicate event name: ''%s''', name);
            end
        end
        reg_ev(nef+1) = e;
    end
end
