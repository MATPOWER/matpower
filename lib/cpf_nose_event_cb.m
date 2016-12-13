function [nx, cx, done, rollback, evnts, cb_data, results] = cpf_nose_event_cb(...
        k, nx, cx, px, done, rollback, evnts, cb_data, cb_args, results)
%CPF_NOSE_EVENT_CB  Callback to handle NOSE events
%   [NX, CX, DONE, ROLLBACK, EVNTS, CB_DATA, RESULTS] = 
%       CPF_NOSE_EVENT_CB(K, NX, CX, PX, DONE, ROLLBACK, EVNTS, ...
%                               CB_DATA, CB_ARGS, RESULTS)
%
%   Callback to handle NOSE events, triggered by event function
%   CPF_NOSE_EVENT to indicate the nose point of the continuation curve.
%
%   This function sets the msg field of the event when the nose point has
%   been found, raises the DONE.flag and sets the DONE.msg.
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

%% initialize return value
stop_at = cb_data.mpopt.cpf.stop_at;

%% handle event
if ~rollback || nx.step == 0
    for i = 1:length(evnts)
        if strcmp(evnts(i).name, 'NOSE') && evnts(i).zero
            if nx.step == 0
                evnts(i).msg = ...
                    sprintf('Nose point eliminated by limit induced bifurcation at %d continuation steps, lambda = %.4g.', k, nx.lam);
            else
                evnts(i).msg = ...
                    sprintf('Reached steady state loading limit in %d continuation steps, lambda = %.4g.', k, nx.lam);
            end

            %% the following conditional is only necessary if we also allow
            %% finding the location of the nose-point without terminating
            if ischar(stop_at) && strcmp(stop_at, 'NOSE');
                done.flag = 1;
                done.msg = evnts(i).msg;
            end
            break;
        end
    end
end
