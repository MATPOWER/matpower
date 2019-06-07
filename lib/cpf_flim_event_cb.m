function [nx, cx, done, rollback, evnts, cb_data, results] = cpf_flim_event_cb(...
        k, nx, cx, px, done, rollback, evnts, cb_data, cb_args, results)
%CPF_FLIM_EVENT_CB  Callback to handle FLIM events
%   [NX, CX, DONE, ROLLBACK, EVNTS, CB_DATA, RESULTS] =
%       CPF_NOSE_EVENT_CB(K, NX, CX, PX, DONE, ROLLBACK, EVNTS, ...
%                               CB_DATA, CB_ARGS, RESULTS)
%
%   Callback to handle FLIM (branch flow limit violation) events,
%   triggered by event function CPF_FLIM_EVENT to indicate the point at which
%   a branch flow limit is reached.
%
%   All branch flows are expected to be within limits for the base case,
%   otherwise the continuation terminates.
%
%   This function sets the msg field of the event when the flow in any branch
%   reaches its limit, raises the DONE.flag and sets the DONE.msg.
%
%   See CPF_DEFAULT_CALLBACK for details of the input and output arguments.

%   MATPOWER
%   Copyright (c) 2016-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Ahmad Abubakar Sadiq, Federal University of Technology Minna, Nigeria
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% base case data
mpc = cb_data.mpc_base;
nb = size(mpc.bus, 1);
i2e_bus = mpc.order.bus.i2e;
f = mpc.branch(:, F_BUS);
t = mpc.branch(:, T_BUS);
SrateA = mpc.branch(:,RATE_A);

%% initialize
if k == 0
    %% compute MVA flows
    Sf = sqrt((mpc.branch(:, PF)).^2 + (mpc.branch(:, QF)).^2);
    St = sqrt((mpc.branch(:, PT)).^2 + (mpc.branch(:, QT)).^2);

    %% check initial lines power flows exceed RATE_A limits
    if any(max(Sf,St) > SrateA)
        done.flag = 1;  %% prepare to terminate
        
        %% find the lines and which lim(s)
        iL = find(max(Sf,St) > SrateA);
        for j = 1:length(iL)
            L = iL(j);
            msg = sprintf('branch flow limit violated in base case: branch %d -- %d exceeds limit of %g MVA\n',...
                i2e_bus(f(L)), i2e_bus(t(L)), SrateA(L));
            
            %% prepare to terminate
            done.flag = 1;
            done.msg = msg;
        end
    end
end

%% skip if finalize or done
if k < 0 || done.flag
    return;
end

%% handle event
for i = 1:length(evnts)
    if strcmp(evnts(i).name, 'FLIM') && evnts(i).zero
        if cb_data.mpopt.verbose > 3
            msg = sprintf('%s\n    ', evnts(i).msg);
        else
            msg = '';
        end
        
        %% find the lines and which lim(s)
        iL = evnts(i).idx;
        for j = 1:length(iL)
            L = iL(j);      %% index of critical branch event of interest
            msg = sprintf('%sbranch flow limit reached\nbranch %d -- %d at limit of %g MVA @ lambda = %.4g, in %d continuation steps',...
                msg, i2e_bus(f(L)), i2e_bus(t(L)), SrateA(L), nx.lam, k);
        end
        
        %% prepare to terminate
        done.flag = 1;
        done.msg = msg;
    end
end