function [cb_state, results] = ...
    cpf_default_callback(k, V_c, lam_c, V_p, lam_p, cb_data, cb_state, cb_args, results)
%CPF_DEFAULT_CALLBACK   Default callback function for CPF
%   [CB_STATE, RESULTS] = ...
%       CPF_DEFAULT_CALLBACK(K, V_C, LAM_C, V_P, LAM_P, ...
%                            CB_DATA, CB_STATE, CB_ARGS, RESULTS)
%
%   Default callback function used by RUNCPF. Takes input from current
%   iteration, returns a user defined state struct and on the final call
%   a results struct.
%
%   Inputs:
%       K - continuation step iteration count
%       V_C - vector of complex bus voltages after K-th corrector step
%       LAM_C - value of LAMBDA after K-th corrector step
%       V_P - vector of complex bus voltages after K-th predictor step
%       LAM_P - value of LAMBDA after K-th predictor step
%       CB_DATA - struct containing potentially useful static data,
%           with the following fields (all based on internal indexing):
%           mpc_base - MATPOWER case struct of base state
%           mpc_target - MATPOWER case struct of target state
%           Sxfr - nb x 1 vector of scheduled transfers in p.u.
%           Ybus - bus admittance matrix
%           Yf - branch admittance matrix, "from" end of branches
%           Yt - branch admittance matrix, "to" end of branches
%           pv - vector of indices of PV buses
%           pq - vector of indices of PQ buses
%           ref - vector of indices of REF buses
%           mpopt - MATPOWER options struct
%       CB_STATE - struct to which the user may add fields containing
%           any information the callback function would like to
%           pass from one invokation to the next (avoiding the
%           following 5 field names which are already used by
%           the default callback: V_p, lam_p, V_c, lam_c, iterations)
%       CB_ARGS - struct specified in MPOPT.cpf.user_callback_args
%       RESULTS - initial value of output struct to be assigned to
%           CPF field of results struct returned by RUNCPF
%
%   Outputs:
%       CB_STATE - updated version of CB_STATE input arg
%       RESULTS - updated version of RESULTS input arg
%
%   This function is called in three different contexts, distinguished
%   as follows:
%   (1) initial - called without RESULTS output arg, with K = 0,
%           after base power flow, before 1st CPF step.
%   (2) iterations - called without RESULTS output arg, with K > 0
%           at each iteration, after predictor-corrector step
%   (3) final - called with RESULTS output arg, after exiting
%           predictor-corrector loop, inputs identical to
%           last iteration call
%
%   See also RUNCPF.

%   MATPOWER
%   Copyright (c) 2013-2015 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   $Id$
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% initialize plotting options
plot_level  = cb_data.mpopt.cpf.plot.level;
plot_bus    = cb_data.mpopt.cpf.plot.bus;
if plot_level
    if isempty(plot_bus)    %% no bus specified
        %% pick PQ bus with largest transfer
        [junk, idx] = max(cb_data.Sxfr(cb_data.pq));
        if isempty(idx) %% or bus 1 if there are none
            idx = 1;
        else
            idx = cb_data.pq(idx(1));
        end
        idx_e = cb_data.mpc_target.order.bus.i2e(idx);
    else
        idx_e = plot_bus;   %% external bus number
        idx = full(cb_data.mpc_target.order.bus.e2i(idx_e));
        if idx == 0
            error('cpf_default_callback: %d is not a valid bus number for MPOPT.cpf.plot.bus', idx_e);
        end
    end
end

%%-----  FINAL call  -----
if nargout == 2
    %% assemble results struct
    results = cb_state;     %% initialize results with final state
    results.max_lam = max(cb_state.lam_c);
    results.iterations = k;

    %% finish final lambda-V nose curve plot
    if plot_level
        %% plot the final nose curve
        plot(cb_state.lam_c', ...
            abs(cb_state.V_c(idx,:))', ...
            '-', 'Color', [0.25 0.25 1]);
        axis([0 max([1;max(cb_state.lam_p);max(cb_state.lam_c)])*1.05 ...
            0 max([1;max(abs(cb_state.V_p(idx)));max(abs(cb_state.V_c(idx)))*1.05])]);
        hold off;
    end
%%-----  INITIAL call  -----
elseif k == 0
    %% initialize state
    cb_state = struct(  'V_p', V_p, ...
                        'lam_p', lam_p, ...
                        'V_c', V_c, ...
                        'lam_c', lam_c, ...
                        'iterations', 0);
    
    %% initialize lambda-V nose curve plot
    if plot_level
        plot(cb_state.lam_p(1), abs(cb_state.V_p(idx,1)), '-', 'Color', [0.25 0.25 1]);
        title(sprintf('Voltage at Bus %d', idx_e));
        xlabel('\lambda');
        ylabel('Voltage Magnitude');
        axis([0 max([1;max(cb_state.lam_p);max(cb_state.lam_c)])*1.05 ...
            0 max([1;max(abs(cb_state.V_p(idx)));max(abs(cb_state.V_c(idx)))*1.05])]);
        hold on;
    end
%%-----  ITERATION call  -----
else
    %% update state
    cb_state.V_p   = [cb_state.V_p V_p];
    cb_state.lam_p = [cb_state.lam_p lam_p];
    cb_state.V_c   = [cb_state.V_c V_c];
    cb_state.lam_c = [cb_state.lam_c lam_c];
    cb_state.iterations    = k;

    %% plot single step of the lambda-V nose curve
    if plot_level > 1
        plot([cb_state.lam_c(k); cb_state.lam_p(k+1)], ...
            [abs(cb_state.V_c(idx,k)); abs(cb_state.V_p(idx,k+1))], ...
            '-', 'Color', 0.85*[1 0.75 0.75]);
        plot([cb_state.lam_p(k+1); cb_state.lam_c(k+1)], ...
            [abs(cb_state.V_p(idx,k+1)); abs(cb_state.V_c(idx,k+1))], ...
            '-', 'Color', 0.85*[0.75 1 0.75]);
        plot(cb_state.lam_p(k+1), abs(cb_state.V_p(idx,k+1)), 'x', ...
            'Color', 0.85*[1 0.75 0.75]);
        plot(cb_state.lam_c(k+1)', ...
            abs(cb_state.V_c(idx,k+1))', ...
            '-o', 'Color', [0.25 0.25 1]);
        axis([0 max([1;max(cb_state.lam_p);max(cb_state.lam_c)])*1.05 ...
            0 max([1;max(abs(cb_state.V_p(idx)));max(abs(cb_state.V_c(idx)))*1.05])]);
        drawnow;
        if plot_level > 2
            pause;
        end
    end
end
