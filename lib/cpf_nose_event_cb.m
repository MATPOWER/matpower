function [cb_state, nn, cc, cb_data, terminate, results] = cpf_nose_event_cb(...
        cont_steps, nn, cc, pp, rollback, critical, terminate, ...
        cb_data, cb_state, cb_args, results)
%CPF_NOSE_EVENT_CB  Event handler for NOSE events
%
%   [CB_STATE, NN, CC, CB_DATA, TERMINATE] = CPF_NOSE_EVENT_CB(CONT_STEPS, ...
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

%% initialize return value
stop_at = cb_data.mpopt.cpf.stop_at;

%% handle event
if ~rollback || nn.step == 0
    for i = 1:length(critical)
        if strcmp(critical(i).name, 'NOSE') && strcmp(critical(i).status, 'ZERO')
            if cb_data.mpopt.verbose
                if nn.step == 0
                    fprintf('\nLimit induced bifurcation at %d continuation steps eliminated nose point\n', cont_steps);
                else
                    fprintf('\nReached steady state loading limit in %d continuation steps\n', cont_steps);
                end
            end

            %% log the event

            %% the following conditional is only necessary if we also allow
            %% finding the location of the nose-point without terminating
            if ischar(stop_at) && strcmp(stop_at, 'NOSE');
                terminate = 1;
            end
            break;
        end
    end
end
