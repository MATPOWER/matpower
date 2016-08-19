function [nn, cc, cb_data, terminate, results] = cpf_target_lam_event_cb(...
        cont_steps, nn, cc, pp, rollback, critical, terminate, ...
        cb_data, cb_args, results)
%CPF_TARGET_LAM_EVENT_CB  Event handler for TARGET_LAM events
%
%   [CB_STATE, NN, CC, CB_DATA, TERMINATE] = CPF_TARGET_LAM_EVENT_CB(CONT_STEPS, ...
%           NN, CC, ROLLBACK, CRITICAL, CB_DATA, CB_STATE, CB_ARGS)
%   
%   Inputs:
%       CONT_STEPS : ...
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

%% skip if initialize, finalize or terminate
if cont_steps <= 0 || terminate
    return;
end

%% make stop_at numerical, 0 = FULL, -Inf = NOSE
stop_at = cb_data.mpopt.cpf.stop_at;
if ischar(stop_at)
    if strcmp(stop_at, 'FULL')
        stop_at = 0;
    else                %% NOSE
        stop_at = -Inf;
    end
end

%% handle event
event_detected = 0;
for i = 1:length(critical)
    if strcmp(critical(i).name, 'TARGET_LAM')
        switch critical(i).status
            case 'ZERO'         %% prepare to terminate
                if cb_data.mpopt.verbose
                    if stop_at == 0     %% FULL
                        fprintf('\nTraced full continuation curve in %d continuation steps\n', cont_steps);
                    else                %% target lambda value
                        fprintf('\nReached desired lambda %g in %d continuation steps\n', ...
                            stop_at, cont_steps);
                    end
                end

                %% log the event

                event_detected = 1;
                terminate = 1;
                break;
            case 'INTERVAL'     %% set step-size & parameterization to terminate next time

                %% log the event

                if stop_at == 0     %% FULL
                    cc.this_step = cc.lam;
                else                %% target lambda value
                    cc.this_step = stop_at - cc.lam;
                end
                cc.this_parm = 1;   %% change to natural parameterization
                event_detected = 1;
                break;
        end
    end
end

%% otherwise, check if predicted lambda of next step will overshoot
%% (by more than remaining distance, to play it safe)
if ~event_detected && ~rollback
    if isempty(nn.this_step)
        step = nn.default_step;
    else
        step = nn.this_step;
    end
    [V0, lam0] = cpf_predictor(nn.V, nn.lam, nn.z, step, cb_data.pv, cb_data.pq);
    if stop_at == 0         %% FULL
        if lam0 < -nn.lam
            nn.this_step = nn.lam;
            nn.this_parm = 1;       %% change to natural parameterization
        end
    elseif stop_at > 0      %% target lambda value
        if lam0 > stop_at + (stop_at - nn.lam)
            nn.this_step = stop_at - nn.lam;
            nn.this_parm = 1;       %% change to natural parameterization
        end
    end
end
