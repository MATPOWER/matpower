function [nx, cx, done, rollback, evnts, cb_data, results] = cpf_target_lam_event_cb(...
        k, nx, cx, px, done, rollback, evnts, cb_data, cb_args, results)
%CPF_TARGET_LAM_EVENT_CB  Callback to handle TARGET_LAM events
%   [NX, CX, DONE, ROLLBACK, EVNTS, CB_DATA, RESULTS] = 
%       CPF_TARGET_LAM_EVENT_CB(K, NX, CX, PX, DONE, ROLLBACK, EVNTS, ...
%                               CB_DATA, CB_ARGS, RESULTS)
%
%   Callback to handle TARGET_LAM events, triggered by event function
%   CPF_TARGET_LAM_EVENT to indicate that a target lambda value has been
%   reached or that the full continuation curve has been traced.
%
%   This function sets the msg field of the event when the target lambda has
%   been found, raises the DONE.flag and sets the DONE.msg. If the current
%   or predicted next step overshoot the target lambda, it adjusts the step
%   size to be exactly what is needed to reach the target, and sets the
%   parameterization for that step to be the natural parameterization.
%
%   See CPF_DEFAULT_CALLBACK for details of the input and output arguments.

%   MATPOWER
%   Copyright (c) 2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% skip if initialize, finalize or done
if k <= 0 || done.flag
    return;
end

%% make stop_at numerical, 0 = FULL, -Inf = NOSE
stop_at = cb_data.mpopt.cpf.stop_at;
verbose = cb_data.mpopt.verbose;
if ischar(stop_at)
    if strcmp(stop_at, 'FULL')
        stop_at = 0;
    else                %% NOSE
        stop_at = -Inf;
    end
end

%% handle event
event_detected = 0;
for i = 1:length(evnts)
    if strcmp(evnts(i).name, 'TARGET_LAM')
        event_detected = 1;
        if evnts(i).zero     %% prepare to terminate
            done.flag = 1;
            if stop_at == 0     %% FULL
                done.msg = sprintf('Traced full continuation curve in %d continuation steps', k);
            else                %% target lambda value
                done.msg = sprintf('Reached desired lambda %g in %d continuation steps', ...
                    stop_at, k);
            end
        else                    %% set step-size & parameterization to terminate next time
            if stop_at == 0     %% FULL
                evnts(i).msg = sprintf('%s\n  step %d expected to overshoot full trace, reduce step size and set natural param', evnts(i).msg, k);
            else                %% target lambda value
                evnts(i).msg = sprintf('%s\n  step %d expected to overshoot target lambda, reduce step size and set natural param', evnts(i).msg, k);
            end
            evnts(i).log = 1;
            if stop_at == 0     %% FULL
                cx.this_step = cx.lam;
            else                %% target lambda value
                cx.this_step = stop_at - cx.lam;
            end
            cx.this_parm = 1;   %% change to natural parameterization
        end
        break;
    end
end

%% otherwise, check if predicted lambda of next step will overshoot
%% (by more than remaining distance, to play it safe)
if ~event_detected && ~rollback
    if isempty(nx.this_step)
        step = nx.default_step;
    else
        step = nx.this_step;
    end
    [V_hat, lam_hat] = cpf_predictor(nx.V, nx.lam, nx.z, step, cb_data.pv, cb_data.pq);
    if stop_at == 0         %% FULL
        if lam_hat < -nx.lam
            nx.this_step = nx.lam;
            nx.this_parm = 1;       %% change to natural parameterization
            if verbose > 2
                fprintf('  step %d expected to overshoot full trace, reduce step size and set natural param\n', k+1);
            end
        end
    elseif stop_at > 0      %% target lambda value
        if lam_hat > stop_at + (stop_at - nx.lam)
            nx.this_step = stop_at - nx.lam;
            nx.this_parm = 1;       %% change to natural parameterization
            if verbose > 2
                fprintf('  step %d expected to overshoot target lambda, reduce step size and set natural param\n', k+1);
            end
        end
    end
end
