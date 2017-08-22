function [nx, cx, done, rollback, evnts, cb_data, results] = cpf_plim_event_cb(...
        k, nx, cx, px, done, rollback, evnts, cb_data, cb_args, results)
%CPF_PLIM_EVENT_CB  Callback to handle PLIM events
%   [NX, CX, DONE, ROLLBACK, EVNTS, CB_DATA, RESULTS] =
%       CPF_PLIM_EVENT_CB(K, NX, CX, PX, DONE, ROLLBACK, EVNTS, ...
%                               CB_DATA, CB_ARGS, RESULTS)
%
%   Callback to handle PLIM (generator active power limit violation) events,
%   triggered by event function CPF_PLIM_EVENT to indicate the point at which
%   an upper active power output limit is reached for a generator.
%
%   When an active power limit is encountered, this function zeros out
%   subsequent transfers from that generator, chooses a new reference bus
%   if necessary, and updates the CB_DATA accordingly, setting the next
%   step size to zero. The event msg is updated with the details of the
%   changes. It also requests termination if all generators reach PMAX.
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

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

mpc = [];

%% handle event
for i = 1:length(evnts)
    if strcmp(evnts(i).name, 'PLIM') && evnts(i).zero
        %% get updated MPC, if necessary
        if isempty(mpc)
            d = cb_data;
            if length(d.ref) ~= 1
                error('cpf_plim_event_cb: ''cpf.enforce_plims'' option only valid for systems with exactly one REF bus')
            end
            mpc = cpf_current_mpc(d.mpc_base, d.mpc_target, ...
                d.Ybus, d.Yf, d.Yt, d.ref, d.pv, d.pq, nx.V, nx.lam, d.mpopt);
            ng = size(mpc.gen, 1);
            i2e_bus = mpc.order.bus.i2e;
            i2e_gen = mpc.order.gen.i2e;
        end

        %% find the generator(s) and which lim(s)
        if cb_data.mpopt.verbose > 3
            msg = sprintf('%s\n    ', evnts(i).msg);
        else
            msg = '';
        end
        ig = evnts(i).idx;
        for j = 1:length(ig)
            g = ig(j);                  %% index of gen of interest
            ib = mpc.gen(g, GEN_BUS);   %% corresponding bus index
            msg = sprintf('%sgen %d @ bus %d reached %g MW Pmax lim @ lambda = %.4g', ...
                msg, i2e_gen(g), i2e_bus(ib), mpc.gen(g, PMAX), nx.lam);
            if ib == cb_data.ref    %% if it is at the ref bus
                %% find gens that are on and at PMAX
                idx_pmax = find( mpc.gen(:, GEN_STATUS) > 0 & ...
                    abs(mpc.gen(:, PG) - mpc.gen(:, PMAX)) < d.mpopt.cpf.p_lims_tol );

                %% find a new ref bus
                candidates = zeros(size(mpc.bus, 1), 1);
                candidates(cb_data.pv) = 1;
                candidates(mpc.gen(idx_pmax, GEN_BUS)) = 0;
                candidates(ib) = 0;
                new_ref = find(candidates, 1);
                if isempty(new_ref)
                    done.flag = 1;
                    done.msg = 'All generators at Pmax';
                else
                    %% convert this bus to PV, set new REF bus
                    mpc.bus(ib, BUS_TYPE) = PV;
                    mpc.bus(new_ref, BUS_TYPE) = REF;

                    %% update new bus types
                    [ref, pv, pq] = bustypes(mpc.bus, mpc.gen);
                    msg = sprintf('%s : ref changed from bus %d to %d', ...
                        msg, i2e_bus(ib), i2e_bus(new_ref));
                end
            end

            %% set Pg to exact limit
            mpc.gen(g, PG) = mpc.gen(g, PMAX);

            %% update new bus types, including in base and target cases
            if ib == cb_data.ref && ~isempty(new_ref)   %% if ref bus has changed
                cb_data.ref = ref;
                cb_data.pv  = pv;
                cb_data.pq  = pq;
                cb_data.mpc_base.bus(  ib, BUS_TYPE) = mpc.bus(ib, BUS_TYPE);
                cb_data.mpc_target.bus(ib, BUS_TYPE) = mpc.bus(ib, BUS_TYPE);
                cb_data.mpc_base.bus(  new_ref, BUS_TYPE) = mpc.bus(new_ref, BUS_TYPE);
                cb_data.mpc_target.bus(new_ref, BUS_TYPE) = mpc.bus(new_ref, BUS_TYPE);
            end

            %% update PG for P limited gen in base and target
            %% (no more active transfer for this gen)
            cb_data.mpc_base.gen(  g, PG) = mpc.gen(g, PG);
            cb_data.mpc_target.gen(g, PG) = mpc.gen(g, PG);
            cb_data.idx_pmax = [cb_data.idx_pmax; g];

            %% update functions
            b = cb_data.mpc_base;
            t = cb_data.mpc_target;
            cb_data.Sbusb = @(Vm)makeSbus(b.baseMVA, b.bus, b.gen, d.mpopt, Vm);
            cb_data.Sbust = @(Vm)makeSbus(t.baseMVA, t.bus, t.gen, d.mpopt, Vm);
            
            %% set size of next step to zero
            nx.this_step = 0;
        end
        evnts(i).msg = msg;
    end
end
