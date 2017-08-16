function [nx, cx, done, rollback, evnts, cb_data, results] = cpf_vlim_event_cb(...
        k, nx, cx, px, done, rollback, evnts, cb_data, cb_args, results)
%CPF_VLIM_EVENT_CB  Callback to handle Voltage limit violation events
%   [NX, CX, DONE, ROLLBACK, EVNTS, CB_DATA, RESULTS] = 
%       CPF_NOSE_EVENT_CB(K, NX, CX, PX, DONE, ROLLBACK, EVNTS, ...
%                               CB_DATA, CB_ARGS, RESULTS)
%
%   Callback to handle VOLTAGE limit events, triggered by event function
%   CPF_VLIM_EVENT to indicate that a bus has reached voltage limit during 
%   the continuation curve. 
%%   It is expected that at base case (mpcb), all initial bus voltages are within limits
%%   (Vmin <= VM <= Vmax), otherwise the continuation terminates. 

%   This function sets the msg field of the event when any of the bus voltages reaches 
%    a limit, raises the DONE.flag and sets the DONE.msg.

%   MATPOWER
%   Copyright (c) 2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Ahmad Abubakar Sadiq, Federal University of Technology Minna, Nigeria
%   and Shrirang Abhyankar, Argonne National Laboratory

%%   This file is not yet part of MATPOWER.
%   It is not yet covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.


%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

%% Initialize
if k <= 0
    d = cb_data;
    %% check initial bus Vm
    if any(d.mpc_base.bus(:, VM) < d.mpc_base.bus(:, VMIN)) || any(d.mpc_base.bus(:, VM) > d.mpc_base.bus(:, VMAX))
        
        mpc = cpf_current_mpc(d.mpc_base, d.mpc_target, ...
            d.Ybus, d.Yf, d.Yt, d.ref, d.pv, d.pq, nx.V, nx.lam, d.mpopt);
        nb = size(mpc.bus, 1);
        i2e_bus = cb_data.mpc_target.order.bus.i2e;

        %% find the bus(es) and which lim(s)
        ib = [find(d.mpc_base.bus(:, VM) < d.mpc_base.bus(:, VMIN));d.mpc_base.bus(:, VM) > d.mpc_base.bus(:, VMAX)];
        for j = 1:length(ib)
            b = ib(j);                 %%index of critical bus event of interest
            minlim = 1;
            if b > nb
                b = b - nb;
                minlim = 0;
            end
            if minlim
                msg = sprintf('Voltage magnitude limits enforced: bus %d exceeds VMIN limit %3.2f at base case',...
                   i2e_bus(b), mpc.bus(b,VMIN));
            else
                msg = sprintf('Voltage magnitude limits enforced: bus %d exceeds VMAX limit %3.2f at base case',...
                   i2e_bus(b), mpc.bus(b,VMAX));
            end
        end 

        done.flag = 1;
        done.msg = msg;   
    end
end
    
%% Finalize or done
if done.flag
    return;
end


%% handle event 
for i = 1:length(evnts)
    if strcmp(evnts(i).name, 'VLIM') && evnts(i).zero
        d = cb_data;
        mpc = cpf_current_mpc(d.mpc_base, d.mpc_target, ...
            d.Ybus, d.Yf, d.Yt, d.ref, d.pv, d.pq, nx.V, nx.lam, d.mpopt);
        nb = size(mpc.bus, 1);
        i2e_bus = cb_data.mpc_target.order.bus.i2e;
        
        %% find the bus(es) and which lim(s)
        if cb_data.mpopt.verbose > 3
            msg = sprintf('%s\n    ', evnts(i).msg);
        else
            msg = '';
        end
        ib = evnts(i).idx;
        for j = 1:length(ib)
            b = ib(j);                 %%index of critical bus event of interest
            minlim = 1;
            if b > nb
                b = b - nb;
                minlim = 0;
            end
            if minlim
                msg = sprintf('%s Bus voltage magnitude limit reached\nBus %d at VMIN lim %3.2f @ lambda = %.4g, in %d continuation steps',...
                       msg, i2e_bus(b), mpc.bus(b,VMIN),nx.lam, k);  %%  i2e_bus(ibb)
            else
                msg = sprintf('%s Bus voltage magnitude limit reached\nbus %d at VMAX lim %3.2f @ lambda = %.4g',...
                       msg, i2e_bus(b), mpc.bus(b,VMAX), nx.lam);  %%  i2e_bus(ibb)
            end
        end 
        
       done.flag = 1;
       done.msg = msg; 
    end
end
