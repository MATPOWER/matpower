function [rollback, ev, nefv] = pne_detect_events(reg_ev, nefv, cefv, step)
% pne_detect_events - Detect events from event function values.
% ::
%
%   [ROLLBACK, CRITICAL_EVENTS, NEFV] = ...
%       PNE_DETECT_EVENTS(REG_EV, NEFV, CEFV, STEP)
%
%   Determines whether any events have been detected based on the
%   corresponding "current" and "next" event function values. If any events
%   are detected, the relevant information is returned in CRITICAL_EVENTS,
%   along with a ROLLBACK flag indicating whether to roll back the proposed
%   next step.
%
%   Inputs:
%       REG_EV : struct containing info about registered event fcns
%       NEFV : cell array of Next Event Function Values
%       CEFV : cell array of Current Event Function Values
%       STEP : current step size
%
%   Outputs:
%       ROLLBACK : flag indicating whether any event has requested a
%           rollback step, if ROLLBACK is true, then CRITICAL_EVENTS will
%           have length 1
%       CRITICAL_EVENTS : struct array containing information about any
%           detected events, with fields:
%           eidx        : event index, index in REG_EV, 0 if no event detected
%           name        : name of event function, empty if none detected
%           zero        : 1 if zero has been detected, 0 otherwise
%                           (i.e. interval detected or no event detected)
%           idx         : index(es) of critical elements in event function
%           step_scale  : linearly interpolated estimate of scaling factor
%                         for current step size required to reach event zero
%           log         : 1 log the event in the results, 0 don't log the event
%                         (set to 1 for zero events, 0 otherwise, can be
%                           modified by callbacks)
%           msg         : event message, set to something generic like
%                           'ZERO detected for TARGET_LAM event' or
%                           'INTERVAL detected for QLIM(3) event', but intended
%                         to be changed/updated by callbacks
%       NEFV : cell array of Next Event Function Values
%           (possibly updated to be exactly zero to avoid re-detection of
%           an interval following a zero detection)

%   MP-Opt-Model
%   Copyright (c) 2016-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% initialize result variables
rollback = 0;
ev = struct( ...            %% initialize struct for storing critical events
    'eidx', 0, ...
    'zero', 0, ...
    'step_scale', 1, ...
    'log', 0, ...
    'name', '', ...
    'idx', 0, ...
    'msg', '' ...
);

%% other initialization
i = 1;              %% index into ev struct
nef = length(nefv); %% number of event functions

%% detect events, first look for event intervals for events requesting rollback
for eidx = 1:nef
    if ~reg_ev(eidx).locate     %% if event doesn't request rollback to locate zero
        continue;               %%   skip to next event
    end

    %% current and previous event function signs
    c_sign = sign(nefv{eidx});
    p_sign = sign(cefv{eidx});

    %% if there's been a sign change and we aren't within event tolerance ...
    idx = find( abs(c_sign) == 1 & c_sign == -p_sign & ...
                abs(nefv{eidx}) > reg_ev(eidx).tol  );
    if ~isempty(idx)
        if step == 0    %% if it's a "repeat" step
            %% (e.g. after warmstart with possible fcn change)
            %% ... make this one the critical one and call it a ZERO event
            ev.eidx = eidx;
            ev.zero = 1;
            ev.step_scale = 1;
            ev.log = 1;
            ev.name = reg_ev(eidx).name;
            ev.idx = idx;
            ev.msg = 'ZERO (BIFURCATION)';
            i = i + 1;
            break;
        else
            %% ... compute step size scaling factors and find index of smallest one
            [step_scale, j] = ...
                min(cefv{eidx}(idx) ./ (cefv{eidx}(idx) - nefv{eidx}(idx)) );

            %% if it's smaller than the current critical one ...
            if step_scale < ev.step_scale
                %% ... make this one the critical one
                ev.eidx = eidx;
                ev.zero = 0;
                ev.step_scale = step_scale;
                ev.log = 0;
                ev.name = reg_ev(eidx).name;
                ev.idx = idx(j);
                ev.msg = 'INTERVAL';
                rollback = 1;   %% signal that a rollback event has been detected
            end
        end
    end
end

%% if no rollback events were detected
if rollback == 0
    %% search for event zeros
    for eidx = 1:nef
        %% if there's an event zero ...
        idx = find( abs(nefv{eidx}) <= reg_ev(eidx).tol );
        if ~isempty(idx)
            %% set event function to exactly zero
            %% (to prevent possible INTERVAL detection again on next step)
            nefv{eidx}(idx) = 0;

            %% ... make this one the critical one
            ev(i).eidx = eidx;
            ev(i).zero = 1;
            ev(i).step_scale = 1;
            ev(i).log = 1;
            ev(i).name = reg_ev(eidx).name;
            ev(i).idx = idx;
            ev(i).msg = 'ZERO';
            i = i + 1;
        end
    end
    
    %% and if no zeros were detected
    if i == 1
        %% search for intervals for non-rollback events
        for eidx = 1:nef
            %% current and previous event function signs
            c_sign = sign(nefv{eidx});
            p_sign = sign(cefv{eidx});

            %% if there's been a sign change ...
            idx = find( abs(c_sign) == 1 & c_sign == -p_sign );
            if ~isempty(idx)
                %% ... compute step size scaling factors ...
                step_scale = cefv{eidx}(idx) ./ (cefv{eidx}(idx) - nefv{eidx}(idx));

                %% ... and save the info as an interval detection
                ev(i).eidx = eidx;
                ev(i).zero = 0;
                ev(i).step_scale = step_scale;
                ev(i).log = 0;
                ev(i).name = reg_ev(eidx).name;
                ev(i).idx = idx;
                ev(i).msg = 'INTERVAL';
                i = i + 1;
            end
        end
    end
end

%% update msgs
if ev(1).eidx
    for i = 1:length(ev)
        if length(nefv{ev(i).eidx}) > 1
            s1 = sprintf('(%d)', ev(i).idx);
        else
            s1 = '';
        end
        if rollback
            s2 = sprintf(' : ROLLBACK by %g', ev(i).step_scale);
        else
            s2 = '';
        end
        ev(i).msg = sprintf('%s detected for %s%s event%s', ...
            ev(i).msg, ev(i).name, s1, s2);
    end
end
