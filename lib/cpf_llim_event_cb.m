function [nx, cx, done, rollback, evnts, cb_data, results] = cpf_llim_event_cb(...
        k, nx, cx, px, done, rollback, evnts, cb_data, cb_args, results)
%CPF_LLIM_EVENT_CB  Callback to handle lines thermal limit violation events
%   [NX, CX, DONE, ROLLBACK, EVNTS, CB_DATA, RESULTS] = 
%       CPF_NOSE_EVENT_CB(K, NX, CX, PX, DONE, ROLLBACK, EVNTS, ...
%                               CB_DATA, CB_ARGS, RESULTS)
%
%   Callback to handle lines thermal limit events, triggered by event function
%   CPF_LLIM_EVENT to indicate that a line has reached its thermal limit during 
%   the continuation curve. 
%%   It is assumed that at base casedata (mpcb), all initial line power flows are within limits
%%   (Sij <= Smax or Sji <= Smax), otherwise an error massage. 

%   This function sets the msg field of the event when any of the lines
%   reach its thermal limits in MVA, raises the DONE.flag and sets the DONE.msg.

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
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% Initialize
if k <= 0
    d = cb_data;
    mpc = cpf_current_mpc(d.mpc_base, d.mpc_target, ...
        d.Ybus, d.Yf, d.Yt, d.ref, d.pv, d.pq, nx.V, nx.lam, d.mpopt);
    f = mpc.branch(:,F_BUS);
    t = mpc.branch(:,T_BUS);
    i2e_bus = cb_data.mpc_target.order.bus.i2e;
    SrateA = d.mpc_base.branch(:,RATE_A);
    Sf = sqrt((mpc.branch(:, PF)).^2 + (mpc.branch(:, QF)).^2);
    St = sqrt((mpc.branch(:, PT)).^2 + (mpc.branch(:, QT)).^2);

    %% check initial lines power flows exceed Rate A limits
    if any(max(Sf,St) > SrateA)
        done.flag = 1; %% Prepare to terminate
        
        %% find the lines and which lim(s)
        iL = find(max(Sf,St) > SrateA);
        for j = 1:length(iL)
            L = iL(j);                          
            msg = sprintf('Line flow limits enforced: Line %d -- %d exceeds MVA limit %3.2f at base case\n',...
                i2e_bus(f(L)), i2e_bus(t(L)), SrateA(L));
               
            %%prepare to terminate
            done.flag = 1;
            done.msg = msg; 
        end
    end
end
    
    
    
if done.flag
    return;
end
%% initialize return value
%stop_at = cb_data.mpopt.cpf.stop_at;
%verbose = cb_data.mpopt.verbose;

%% specify Voltage limit tolerance
%Vlim_tol = d.mpopt.cpf.v_lims_tol;



mpc = [];

%% handle event 
for i = 1:length(evnts)
    if strcmp(evnts(i).name, 'LLIM') && evnts(i).zero
        %% get updated MPC, if necessary
        if isempty (mpc)
            d = cb_data;
            mpc = cpf_current_mpc(d.mpc_base, d.mpc_target, ...
                d.Ybus, d.Yf, d.Yt, d.ref, d.pv, d.pq, nx.V, nx.lam, d.mpopt);
            f = mpc.branch(:,F_BUS);
            t = mpc.branch(:,T_BUS);
            i2e_bus = cb_data.mpc_target.order.bus.i2e;
            SrateA = d.mpc_base.branch(:,RATE_A);
            Sf = sqrt((mpc.branch(:, PF)).^2 + (mpc.branch(:, QF)).^2);
            St = sqrt((mpc.branch(:, PT)).^2 + (mpc.branch(:, QT)).^2);
                                   
            %% check initial lines power flows Sf and St
            if any(max(Sf,St) > SrateA)
                done.flag = 1;
                error('Base case power flow violates line limits');
            end
        end
        %% find the lines and which lim(s)
        if cb_data.mpopt.verbose > 3
            msg = sprintf('%s\n    ', evnts(i).msg);
        else
            msg = '';
        end
        iL = evnts(i).idx;
        for j = 1:length(iL)
            L = iL(j);                                  %%index of critical line event of interest
            msg = sprintf('%s Line flow limit reached\nLine %d -- %d at MVA limit %3.2f @ lambda = %.4g, in %d continuation steps',...
                msg, i2e_bus(f(L)), i2e_bus(t(L)), SrateA(L), nx.lam, k);  %%  i2e_bus(ibb)
               
            %%prepare to terminate
            done.flag = 1;
            done.msg = msg; 
        end
    end
%                else
%                    if (Sij - maxlim) <= d.mpopt.cpf.L_lims_tol             %%prepare to terminate
%                     done.flag = 1;
%                     done.msg = sprintf('%s Line %d from bus %d to bus %d, limits violation within Llim_tol in %d continuation steps, lambda = %.4g.', msg, L, i2e_bus(f(L)), i2e_bus(t(L)), k, nx.lam);
%                    elseif (Sji - maxlim) <= d.mpopt.cpf.v_lims_tol
%                        done.flag = 1;
%                         done.msg = sprintf('%s Line %d from bus %d to bus %d, limits violation within Llim_tol in %d continuation steps, lambda = %.4g.', msg, L, i2e_bus(t(L)), i2e_bus(f(L)), k, nx.lam);
%                    end
%                    %% check if any current line exceed thermal limits in MVA
%                    if Sij > mpc.branch(:, RATE_A)
%                        evnts(i).msg = sprintf('%s\n step expected to violate Llim_tol, reduce step size and set natural param', evnts(i), k);
%                        evnts(i).log = 1;
%                    elseif  Sji > mpc.branch(:, RATE_A)
%                        evnts(i).msg = sprintf('%s\n step expected to violate Llim_tol, reduce step size and set natural param', evnts(i), k);
%                        evnts(i).log = 1;
%                    end
%                    cx.this_step = cx.default_step;
%                    cx.this_parm = 1; %% change to natural parameterization or reduce step size
%                end
%            end
%            %% otherwise, check if the predicted complex power flows Sij or Sji of next step will violate Llim
%            if ~event_detected && ~rollback
%                if isempty(nx.this_step)
%                    step = nx.default_step;
%                else
%                    step = nx.this_step;
%                end
%                [V_hat, lam_hat] = cpf_predictor(nx.V, nx.lam, nx.z, step, cb_data.pv, cb_data.pq);
%                if Sij > Smax || Sji > Smax
%                    nx.this_parm = 1; %% change to natural parameterization
%                    if verbose > 2
%                        fprintf(' step %d expected to violate Vlim, reduce step and set natural param\n', k+1);
%                    end
%                end
%            end
%            break
end