function [rollback, critical, cef] = cpf_detect_events(cpf_events, cef, pef, step, verbose)
%CPF_DETECT_EVENTS  Detect events from event function values
%
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
%           k           : index of event in list of registered events
%           name        : name of event function
%           idx         : index(es) of critical elements in event function
%           step_scale  : linearly interpolated estimate of scaling factor
%                         for current step size required to reach event zero
%           status      : type of event detected 'INTERVAL' or 'ZERO'
%       CEF : cell array of Current Event Function values

%   MATPOWER
%   Copyright (c) 2016 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% initialize result variables
rollback = 0;
critical = struct('k', 0, 'name', '', 'idx', 0, 'step_scale', 1, 'status', '');

%% other initialization
i = 1;              %% index into critical struct
nef = length(cef);  %% number of event functions

%% detect events, first look for event intervals for events requesting rollback
for k = 1:nef
    if ~cpf_events(k).rollback  %% if event does not request rollback
        continue;               %%   skip to next event
    end

    %% current and previous event function signs
    c_sign = sign(cef{k});
    p_sign = sign(pef{k});

    %% if there's been a sign change and we aren't within event tolerance ...
    idx = find( abs(c_sign) == 1 & c_sign == -p_sign & ...
                abs(cef{k}) > cpf_events(k).tol  );
    if ~isempty(idx)
        if step == 0    %% if it's a "repeat" step (e.g. after bus type changes)
            %% ... make this one the critical one and call it a ZERO event
            critical.k = k;
            critical.name = cpf_events(k).name;
            critical.idx = idx;
            critical.step_scale = 1;
            critical.status = 'ZERO';
            i = i + 1;
            break;
        else
            %% ... compute step size scaling factors and find index of smallest one
            [step_scale, j] = ...
                min(pef{k}(idx) ./ (pef{k}(idx) - cef{k}(idx)) );

            %% if it's smaller than the current critical one ...
            if step_scale < critical.step_scale
                %% ... make this one the critical one
                critical.k = k;
                critical.name = cpf_events(k).name;
                critical.idx = idx(j);
                critical.step_scale = step_scale;
                critical.status = 'INTERVAL';
                rollback = 1;   %% signal that a rollback event has been detected
            end
        end
    end
end

%% if no rollback events were detected
if rollback == 0
    %% search for event zeros
    for k = 1:nef
        %% if there's an event zero ...
        idx = find( abs(cef{k}) <= cpf_events(k).tol );
        if ~isempty(idx)
            %% set event function to exactly zero
            %% (to prevent possible INTERVAL detection again on next step)
            cef{k}(idx) = 0;

            %% ... make this one the critical one
            critical(i).k = k;
            critical(i).name = cpf_events(k).name;
            critical(i).idx = idx;
            critical(i).step_scale = 1;
            critical(i).status = 'ZERO';
            i = i + 1;
        end
    end
    
    %% and if no zeros were detected
    if i == 1
        %% search for intervals for non-rollback events
        for k = 1:nef
            %% current and previous event function signs
            c_sign = sign(cef{k});
            p_sign = sign(pef{k});

            %% if there's been a sign change ...
            idx = find( abs(c_sign) == 1 & c_sign == -p_sign );
            if ~isempty(idx)
                %% ... compute step size scaling factors ...
                step_scale = pef{k}(idx) ./ (pef{k}(idx) - cef{k}(idx));

                %% ... and save the info
                critical(i).k = k;
                critical(i).name = cpf_events(k).name;
                critical(i).idx = idx;
                critical(i).step_scale = step_scale;
                critical(i).status = 'INTERVAL';
                i = i + 1;
            end
        end
    end
end
if verbose > 2 && critical(1).k
    for i = 1:length(critical)
        ce = critical(i);
        k = ce.k;
        fprintf('   %s detected for %s event : ', ce.status, ce.name);
        if rollback
            fprintf('ROLLBACK by %g\n', ce.step_scale);
        else
            fprintf('CONTINUE\n');
        end
    end
end
