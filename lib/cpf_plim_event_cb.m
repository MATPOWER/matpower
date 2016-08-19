function [cb_state, nn, cc, cb_data, terminate, results] = cpf_plim_event_cb(...
        cont_steps, nn, cc, pp, rollback, critical, terminate, ...
        cb_data, cb_state, cb_args, results)
%CPF_PLIM_EVENT_CB  Event handler for NOSE events
%
%   [CB_STATE, NN, CC, CB_DATA, TERMINATE] = CPF_PLIM_EVENT_CB(CONT_STEPS, ...
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

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

mpc = [];

%% handle event
for i = 1:length(critical)
    if strcmp(critical(i).name, 'PLIM') && strcmp(critical(i).status, 'ZERO')
        %% get updated MPC, if necessary
        if isempty(mpc)
            d = cb_data;
            if length(d.ref) ~= 1
                error('cpf_plim_event_cb: ''cpf.enforce_plims'' option only valid for systems with exactly one REF bus')
            end
            mpc = cpf_current_mpc(d.mpc_base, d.mpc_target, ...
                d.Ybus, d.Yf, d.Yt, d.ref, d.pv, d.pq, nn.V, nn.lam, d.mpopt);
            ng = size(mpc.gen, 1);
            i2e_bus = cb_data.mpc_target.order.bus.i2e;
            i2e_gen = cb_data.mpc_target.order.gen.i2e;
        end
        
        %% find the generator(s) and which lim(s)
        ig = critical(i).idx;
        for k = 1:length(ig)    %% ig(k) is the gen index of interest
            ib = mpc.gen(ig(k), GEN_BUS);   %% corresponding bus index
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
                    terminate = 1;
                else
                    %% convert this bus to PV, set new REF bus
                    mpc.bus(ib, BUS_TYPE) = PV;
                    mpc.bus(new_ref, BUS_TYPE) = REF;

                    %% update new bus types
                    [ref, pv, pq] = bustypes(mpc.bus, mpc.gen);
                end
            end

            if cb_data.mpopt.verbose
                if ib == cb_data.ref   %% at ref bus
                    if isempty(new_ref)
                        fprintf('  gen %d @ bus %d reached %g MW Pmax lim @ lambda = %.3g\nCPF Termination : all generators at Pmax\n', ...
                            i2e_gen(ig(k)), i2e_bus(ib), mpc.gen(ig(k), PMAX), nn.lam);
                    else
                        fprintf('  gen %d @ bus %d reached %g MW Pmax lim @ lambda = %.3g : ref changed from bus %d to %d\n', ...
                            i2e_gen(ig(k)), i2e_bus(ib), mpc.gen(ig(k), PMAX), nn.lam, ...
                            i2e_bus(ib), i2e_bus(new_ref));
                    end
                else
                    fprintf('  gen %d @ bus %d reached %g MW Pmax lim @ lambda = %.3g\n', ...
                        i2e_gen(ig(k)), i2e_bus(ib), mpc.gen(ig(k), PMAX), nn.lam);
                end
            end

            %% log the event
            

            %% set Pg to exact limit
            mpc.gen(ig(k), PG) = mpc.gen(ig(k), PMAX);

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
            cb_data.mpc_base.gen(  ig(k), PG) = mpc.gen(ig(k), PG);
            cb_data.mpc_target.gen(ig(k), PG) = mpc.gen(ig(k), PG);
            cb_data.idx_pmax = [cb_data.idx_pmax; ig(k)];

            %% update functions
            b = cb_data.mpc_base;
            t = cb_data.mpc_target;
            cb_data.Sbusb = @(Vm)makeSbus(b.baseMVA, b.bus, b.gen, d.mpopt, Vm);
            cb_data.Sbust = @(Vm)makeSbus(t.baseMVA, t.bus, t.gen, d.mpopt, Vm);
            cb_data.Sxfr  = @(Vm)(cb_data.Sbust(Vm) - cb_data.Sbusb(Vm));
            
            %% set size of next step to zero
            nn.this_step = 0;
        end
    end
end
