function [cb_state, nn, cc, cb_data, terminate, results] = cpf_qlim_event_cb(...
        cont_steps, nn, cc, pp, rollback, critical, terminate, ...
        cb_data, cb_state, cb_args, results)
%CPF_QLIM_EVENT_CB  Event handler for NOSE events
%
%   [CB_STATE, NN, CC, CB_DATA, TERMINATE] = CPF_QLIM_EVENT_CB(CONT_STEPS, ...
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
    if strcmp(critical(i).name, 'QLIM') && strcmp(critical(i).status, 'ZERO')
        %% get updated MPC, if necessary
        if isempty(mpc)
            d = cb_data;
            if length(d.ref) ~= 1
                error('cpf_qlim_event_cb: ''cpf.enforce_qlims'' option only valid for systems with exactly one REF bus')
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
            maxlim = 1;
            if ig(k) > ng
                ig(k) = ig(k) - ng;
                maxlim = 0;
            end

            if cb_data.mpopt.verbose
                if maxlim
                    fprintf('  gen %d @ bus %d reached %g MVAr Qmax lim @ lambda = %.3g : bus %d converted to PQ\n', ...
                        i2e_gen(ig(k)), i2e_bus(ib), mpc.gen(ig(k), QMAX), nn.lam, i2e_bus(ib));
                else
                    fprintf('  gen %d @ bus %d reached %g MVAr Qmin lim @ lambda = %.3g : bus %d converted to PQ\n', ...
                        i2e_gen(ig(k)), i2e_bus(ib), mpc.gen(ig(k), QMIN), nn.lam, i2e_bus(ib));
                end
            end

            %% log the event
            

            %% set Qg to exact limit and convert the generator's bus to PQ bus
            if maxlim
                mpc.gen(ig(k), QG) = mpc.gen(ig(k), QMAX);
            else
                mpc.gen(ig(k), QG) = mpc.gen(ig(k), QMIN);
            end
            mpc.bus(ib, BUS_TYPE) = PQ;

            %% infeasibility check
            on = find(mpc.gen(:, GEN_STATUS) > 0 & ...  %% which generators are on?
                      mpc.bus(mpc.gen(:, GEN_BUS), BUS_TYPE) ~= PQ);  %% ... and are not PQ buses

            if isempty(on)
                terminate = 1;
                if cb_data.mpopt.verbose
                    fprintf('CPF termination: No REF, PV buses remaining\n');
                end
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
                cb_data.mpc_base.gen(ig(k), QG)   = mpc.gen(ig(k), QG);
                cb_data.mpc_target.gen(ig(k), QG) = mpc.gen(ig(k), QG);

                %% update PG for previous slack gen in base and target
                %% (no more active transfer for this gen)
                %% RDZ: Why?
                if oldref ~= ref
                    cb_data.mpc_base.gen(ig(k), PG)   = mpc.gen(ig(k),PG);
                    cb_data.mpc_target.gen(ig(k), PG) = mpc.gen(ig(k),PG);
                end
                
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
end
