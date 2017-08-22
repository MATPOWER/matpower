function [nx, cx, done, rollback, evnts, cb_data, results] = cpf_qlim_event_cb(...
        k, nx, cx, px, done, rollback, evnts, cb_data, cb_args, results)
%CPF_QLIM_EVENT_CB  Callback to handle QLIM events
%   [NX, CX, DONE, ROLLBACK, EVNTS, CB_DATA, RESULTS] =
%       CPF_QLIM_EVENT_CB(K, NX, CX, PX, DONE, ROLLBACK, EVNTS, ...
%                               CB_DATA, CB_ARGS, RESULTS)
%
%   Callback to handle QLIM (generator reactive power limit violation) events,
%   triggered by event function CPF_QLIM_EVENT to indicate the point at which
%   an upper or lower reactive power output limit is reached for a generator.
%
%   When a reactive power limit is encountered, this function zeros out
%   subsequent transfers from that generator, changes it's bus type to PQ,
%   chooses a new reference bus if necessary, and updates the CB_DATA
%   accordingly, setting the next step size to zero. The event msg is updated
%   with the details of the changes. It also requests termination if no more
%   PV or REF buses remain.
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
    if strcmp(evnts(i).name, 'QLIM') && evnts(i).zero
        %% get updated MPC, if necessary
        if isempty(mpc)
            d = cb_data;
            if length(d.ref) ~= 1
                error('cpf_qlim_event_cb: ''cpf.enforce_qlims'' option only valid for systems with exactly one REF bus')
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
            maxlim = 1;
            if g > ng
                g = g - ng;
                maxlim = 0;
            end
            ib = mpc.gen(g, GEN_BUS);   %% corresponding bus index
            if maxlim
                msg = sprintf('%sgen %d @ bus %d reached %g MVAr Qmax lim @ lambda = %.4g : bus %d converted to PQ', ...
                    msg, i2e_gen(g), i2e_bus(ib), mpc.gen(g, QMAX), nx.lam, i2e_bus(ib));
            else
                msg = sprintf('%sgen %d @ bus %d reached %g MVAr Qmin lim @ lambda = %.4g : bus %d converted to PQ', ...
                    msg, i2e_gen(g), i2e_bus(ib), mpc.gen(g, QMIN), nx.lam, i2e_bus(ib));
            end

            %% set Qg to exact limit and convert the generator's bus to PQ bus
            if maxlim
                mpc.gen(g, QG) = mpc.gen(g, QMAX);
            else
                mpc.gen(g, QG) = mpc.gen(g, QMIN);
            end
            mpc.bus(ib, BUS_TYPE) = PQ;

            %% infeasibility check
            on = find(mpc.gen(:, GEN_STATUS) > 0 & ...  %% which generators are on?
                      mpc.bus(mpc.gen(:, GEN_BUS), BUS_TYPE) ~= PQ);  %% ... and are not PQ buses

            if isempty(on)
                done.flag = 1;
                done.msg = 'No REF or PV buses remaining.';
            else
                oldref = cb_data.ref;   %% save previous ref bus
                [ref, pv, pq] = bustypes(mpc.bus, mpc.gen);
                if oldref ~= ref        %% ref bus changed
                    mpc.bus(ref, BUS_TYPE) = REF;
                end

                %% update new bus types, including in base and target cases
                cb_data.ref = ref;
                cb_data.pv  = pv;
                cb_data.pq  = pq;
                cb_data.mpc_base.bus(  :, BUS_TYPE) = mpc.bus(:, BUS_TYPE);
                cb_data.mpc_target.bus(:, BUS_TYPE) = mpc.bus(:, BUS_TYPE);
            
                %% update QG for Q limited gen in base and target
                %% (no more reactive transfer for this gen)
                cb_data.mpc_base.gen(  g, QG) = mpc.gen(g, QG);
                cb_data.mpc_target.gen(g, QG) = mpc.gen(g, QG);

                %% update PG for previous slack gen in base and target
                %% (no more active transfer for this gen)
                if oldref ~= ref
                    cb_data.mpc_base.gen(  g, PG) = mpc.gen(g,PG);
                    cb_data.mpc_target.gen(g, PG) = mpc.gen(g,PG);
                end
                
                %% update functions
                b = cb_data.mpc_base;
                t = cb_data.mpc_target;
                cb_data.Sbusb = @(Vm)makeSbus(b.baseMVA, b.bus, b.gen, d.mpopt, Vm);
                cb_data.Sbust = @(Vm)makeSbus(t.baseMVA, t.bus, t.gen, d.mpopt, Vm);
                
                %% set size of next step to zero
                nx.this_step = 0;
            end
        end
        evnts(i).msg = msg;
    end
end
