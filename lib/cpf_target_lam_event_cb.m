function [nx, cx, done, rollback, evnts, cb_data, results] = cpf_target_lam_event_cb(...
        k, nx, cx, px, done, rollback, evnts, cb_data, cb_args, results)
%CPF_TARGET_LAM_EVENT_CB  Event handler for TARGET_LAM events
%
%   [CB_STATE, NX, CX, CB_DATA, DONE] = CPF_TARGET_LAM_EVENT_CB(CONT_STEPS, ...
%           NX, CX, ROLLBACK, CRITICAL, CB_DATA, CB_STATE, CB_ARGS)
%   
%   Inputs:
%       K : ...
%
%   Outputs:
%       CB_STATE : ...

%   MATPOWER
%   Copyright (c) 2016 by Power System Engineering Research Center (PSERC)
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
