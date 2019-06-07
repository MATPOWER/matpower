function [nx, cx, done, rollback, evnts, cb_data, results] = cpf_vlim_event_cb(...
        k, nx, cx, px, done, rollback, evnts, cb_data, cb_args, results)
%CPF_VLIM_EVENT_CB  Callback to handle VLIM events
%   [NX, CX, DONE, ROLLBACK, EVNTS, CB_DATA, RESULTS] =
%       CPF_VLIM_EVENT_CB(K, NX, CX, PX, DONE, ROLLBACK, EVNTS, ...
%                               CB_DATA, CB_ARGS, RESULTS)
%
%   Callback to handle VLIM (bus voltage magnitude limit violation) events,
%   triggered by event function CPF_VLIM_EVENT to indicate the point at which
%   an upper or lower voltage magnitude limit is reached for a bus.
%
%   All bus voltages are expected to be within limits for the base case,
%   otherwise the continuation terminates.
%
%   This function sets the msg field of the event when the voltage magnitude
%   at any bus reaches its upper or lower limit, raises the DONE.flag and sets
%   the DONE.msg.
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

%% base case data
mpc = cb_data.mpc_base;
nb = size(mpc.bus, 1);
i2e_bus = mpc.order.bus.i2e;

%% initialize
if k == 0
    %% check initial bus Vm
    if any(mpc.bus(:, VM) < mpc.bus(:, VMIN)) || any(mpc.bus(:, VM) > mpc.bus(:, VMAX))
        %% find the bus(es) and which lim(s)
        ib = find([mpc.bus(:, VM) < mpc.bus(:, VMIN);mpc.bus(:, VM) > mpc.bus(:, VMAX)]);
        for j = 1:length(ib)
            b = ib(j);          %% index of critical bus event of interest
            if b > nb
                b = b - nb;
                msg = sprintf('bus voltage magnitude limit violated in base case: bus %d exceeds Vmax limit %g p.u.',...
                   i2e_bus(b), mpc.bus(b,VMAX));
            else
                msg = sprintf('bus voltage magnitude limit violated in base case: bus %d exceeds Vmin limit %g p.u.',...
                   i2e_bus(b), mpc.bus(b,VMIN));
            end
        end
        
        done.flag = 1;
        done.msg = msg;
    end
end

%% skip if finalize or done
if k < 0 || done.flag
    return;
end

%% handle event
for i = 1:length(evnts)
    if strcmp(evnts(i).name, 'VLIM') && evnts(i).zero
        if cb_data.mpopt.verbose > 3
            msg = sprintf('%s\n    ', evnts(i).msg);
        else
            msg = '';
        end
        
        %% find the bus(es) and which lim(s)
        ib = evnts(i).idx;
        for j = 1:length(ib)
            b = ib(j);          %% index of critical bus event of interest
            if b > nb
                b = b - nb;
                msg = sprintf('%sbus voltage magnitude limit reached\nbus %d at VMAX limit %g p.u. @ lambda = %.4g, in %d continuation steps',...
                    msg, i2e_bus(b), mpc.bus(b,VMAX), nx.lam, k);
            else
                msg = sprintf('%sbus voltage magnitude limit reached\nbus %d at Vmin limit %g p.u. @ lambda = %.4g, in %d continuation steps',...
                    msg, i2e_bus(b), mpc.bus(b,VMIN),nx.lam, k);
            end
        end
        
        done.flag = 1;
        done.msg = msg;
    end
end
