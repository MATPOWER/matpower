function [rollback, evnts, cef] = cpf_detect_events(cpf_events, cef, pef, step, verbose)
%CPF_DETECT_EVENTS  Detect events from event function values
%   [ROLLBACK, CRITICAL_EVENTS, CEF] = CPF_DETECT_EVENTS(CPF_EVENTS, CEF, PEF, STEP, VERBOSE)
%   
%   Inputs:
%       CPF_EVENTS : struct containing info about registered CPF event fcns
%       CEF : cell array of Current Event Function values
%       PEF : cell array of Previous Event Function values
%       STEP : current step size
%       VERBOSE : 0 = no output, otherwise level of verbose output
%
%   Outputs:
%       ROLLBACK : flag indicating whether any event has requested a
%           rollback step
%       CRITICAL_EVENTS : struct array containing information about any
%           detected events, with fields:
%           eidx        : event index, in list of registered events
%                           0 if no event detected
%           name        : name of event function, empty if none detected
%           zero        : 1 if zero has been detected, 0 otherwise
%                           (interval detected or no event detected)
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
%       CEF : cell array of Current Event Function values

%   MATPOWER
%   Copyright (c) 2016-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% initialize result variables
rollback = 0;
evnts = struct( ...
    'eidx', 0, ...
    'name', '', ...
    'zero', 0, ...
    'idx', 0, ...
    'step_scale', 1, ...
    'log', 0, ...
    'msg', '' ...
);

%% other initialization
i = 1;              %% index into evnts struct
nef = length(cef);  %% number of event functions

%% detect events, first look for event intervals for events requesting rollback
for eidx = 1:nef
    if ~cpf_events(eidx).locate %% if event doesn't request rollback to locate zero
        continue;               %%   skip to next event
    end

    %% current and previous event function signs
    c_sign = sign(cef{eidx});
    p_sign = sign(pef{eidx});

    %% if there's been a sign change and we aren't within event tolerance ...
    idx = find( abs(c_sign) == 1 & c_sign == -p_sign & ...
                abs(cef{eidx}) > cpf_events(eidx).tol  );
    if ~isempty(idx)
        if step == 0    %% if it's a "repeat" step (e.g. after bus type changes)
            %% ... make this one the critical one and call it a ZERO event
            evnts.eidx = eidx;
            evnts.name = cpf_events(eidx).name;
            evnts.zero = 1;
            evnts.idx = idx;
            evnts.step_scale = 1;
            evnts.log = 1;
            evnts.msg = 'ZERO (BIFURCATION)';
            i = i + 1;
            break;
        else
            %% ... compute step size scaling factors and find index of smallest one
            [step_scale, j] = ...
                min(pef{eidx}(idx) ./ (pef{eidx}(idx) - cef{eidx}(idx)) );

            %% if it's smaller than the current critical one ...
            if step_scale < evnts.step_scale
                %% ... make this one the critical one
                evnts.eidx = eidx;
                evnts.name = cpf_events(eidx).name;
                evnts.zero = 0;
                evnts.idx = idx(j);
                evnts.step_scale = step_scale;
                evnts.log = 0;
                evnts.msg = 'INTERVAL';
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
        idx = find( abs(cef{eidx}) <= cpf_events(eidx).tol );
        if ~isempty(idx)
            %% set event function to exactly zero
            %% (to prevent possible INTERVAL detection again on next step)
            cef{eidx}(idx) = 0;

            %% ... make this one the critical one
            evnts(i).eidx = eidx;
            evnts(i).name = cpf_events(eidx).name;
            evnts(i).zero = 1;
            evnts(i).idx = idx;
            evnts(i).step_scale = 1;
            evnts(i).log = 1;
            evnts(i).msg = 'ZERO';
            i = i + 1;
        end
    end
    
    %% and if no zeros were detected
    if i == 1
        %% search for intervals for non-rollback events
        for eidx = 1:nef
            %% current and previous event function signs
            c_sign = sign(cef{eidx});
            p_sign = sign(pef{eidx});

            %% if there's been a sign change ...
            idx = find( abs(c_sign) == 1 & c_sign == -p_sign );
            if ~isempty(idx)
                %% ... compute step size scaling factors ...
                step_scale = pef{eidx}(idx) ./ (pef{eidx}(idx) - cef{eidx}(idx));

                %% ... and save the info as an interval detection
                evnts(i).eidx = eidx;
                evnts(i).name = cpf_events(eidx).name;
                evnts(i).zero = 0;
                evnts(i).idx = idx;
                evnts(i).step_scale = step_scale;
                evnts(i).log = 0;
                evnts(i).msg = 'INTERVAL';
                i = i + 1;
            end
        end
    end
end

%% update msgs
if evnts(1).eidx
    for i = 1:length(evnts)
        if length(cef{evnts(i).eidx}) > 1
            s1 = sprintf('(%d)', evnts(i).idx);
        else
            s1 = '';
        end
        if rollback
            s2 = sprintf('ROLLBACK by %g', evnts(i).step_scale);
        else
            s2 = 'CONTINUE';
        end
        evnts(i).msg = sprintf('%s detected for %s%s event : %s', ...
            evnts(i).msg, evnts(i).name, s1, s2);
    end
end
