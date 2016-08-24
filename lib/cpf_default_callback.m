function [nx, cx, done, rollback, evnts, cb_data, results] = cpf_default_callback(...
        k, nx, cx, px, done, rollback, evnts, cb_data, cb_args, results)
% function [cb_state, results] = ...
%     cpf_default_callback(k, step, V, lam, V_hat, lam_hat, cb_data, cb_state, cb_args, results)
%CPF_DEFAULT_CALLBACK   Default callback function for CPF
%   [CB_STATE, RESULTS] = ...
%       CPF_DEFAULT_CALLBACK(K, STEP, V, LAM, V_HAT, LAM_HAT, ...
%                            CB_DATA, CB_STATE, CB_ARGS, RESULTS)
%
%   Default callback function used by RUNCPF. Takes input from current
%   iteration, returns a user defined state struct and on the final call
%   a results struct.
%
%   Inputs:
%       K - continuation step iteration count
%       STEP - step size for K-th step
%       V - vector of complex bus voltages after K-th corrector step
%       LAM - value of LAMBDA after K-th corrector step
%       V_HAT - vector of complex bus voltages after K-th predictor step
%       LAM_HAT - value of LAMBDA after K-th predictor step
%       CB_DATA - struct containing potentially useful static data,
%           with the following fields (all based on internal indexing):
%           mpc_base - MATPOWER case struct of base state
%           mpc_target - MATPOWER case struct of target state
%           Sbusb - handle of function returning nb x 1 vector of complex
%               base case injections in p.u. and derivatives w.r.t. |V|
%           Sbust - handle of function returning nb x 1 vector of complex
%               target case injections in p.u. and derivatives w.r.t. |V|
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
%           the default callback: V_hat, lam_hat, V, lam, iterations)
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
%   Copyright (c) 2013-2016 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% skip if rollback, except if it is a FINAL call
if rollback && k > 0
    return;
end

%% initialize variables
step    = nx.step;
V       = nx.V;
lam     = nx.lam;
V_hat   = nx.V_hat;
lam_hat = nx.lam_hat;

%%-----  initialize/update state/results  -----
if k == 0       %% INITIAL call
    %% initialize state
    cxx = struct(   'V_hat',    V_hat, ...
                    'lam_hat',  lam_hat, ...
                    'V',        V, ...
                    'lam',      lam, ...
                    'steps',    step, ...
                    'iterations', 0     );
    nxx = cxx;
    cx.cb.default = cxx;    %% update current callback state
    nx.cb.default = nxx;    %% updatenext callback state
else
    nxx = nx.cb.default;    %% get next callback state
    if k > 0    %% ITERATION call
        %% update state
        nxx.V_hat   = [nxx.V_hat    V_hat];
        nxx.lam_hat = [nxx.lam_hat  lam_hat];
        nxx.V       = [nxx.V        V];
        nxx.lam     = [nxx.lam      lam];
        nxx.steps   = [nxx.steps    step];
        nxx.iterations = k;
        nx.cb.default = nxx;    %% update next callback state
    else            %% FINAL call
        %% assemble results struct
        results.V_hat       = nxx.V_hat;
        results.lam_hat     = nxx.lam_hat;
        results.V           = nxx.V;
        results.lam         = nxx.lam;
        results.steps       = nxx.steps;
        results.iterations  = -k;
        results.max_lam     = max(nxx.lam);
    end
end

%%-----  plot continuation curve  -----
%% initialize plotting options
plot_level  = cb_data.mpopt.cpf.plot.level;
plot_bus    = cb_data.mpopt.cpf.plot.bus;
plot_bus_default = 0;
if plot_level
    if isempty(plot_bus) && ~isfield(nxx, 'plot_bus_default')   %% no bus specified
        %% pick PQ bus with largest transfer
        Sxfr = cb_data.Sbust(abs(V)) - cb_data.Sbusb(abs(V));
        [junk, idx] = max(Sxfr(cb_data.pq));
        if isempty(idx) %% or bus 1 if there are none
            idx = 1;
        else
            idx = cb_data.pq(idx(1));
        end
        idx_e = cb_data.mpc_target.order.bus.i2e(idx);
        
        %% save it to keep it from changing in subsequent calls
        plot_bus_default = idx_e;
    else
        if isempty(plot_bus)
            idx_e = nxx.plot_bus_default;   %% external bus number, saved
        else
            idx_e = plot_bus;               %% external bus number, provided
        end
        idx = full(cb_data.mpc_target.order.bus.e2i(idx_e));
        if idx == 0
            error('cpf_default_callback: %d is not a valid bus number for MPOPT.cpf.plot.bus', idx_e);
        end
    end

    %% set bounds for plot axes
    xmin = 0;
    xmax = max([max(nxx.lam_hat); max(nxx.lam)]);
    ymin = min([min(abs(nxx.V_hat(idx, :))); min(abs(nxx.V(idx, :)))]);
    ymax = max([max(abs(nxx.V_hat(idx, :))); max(abs(nxx.V(idx, :)))]);
    if xmax < xmin + cb_data.mpopt.cpf.step / 100;
        xmax = xmin + cb_data.mpopt.cpf.step / 100;
    end
    if ymax - ymin < 2e-5;
        ymax = ymax + 1e-5;
        ymin = ymin - 1e-5;
    end
    xmax = xmax * 1.05;
    ymax = ymax + 0.05 * (ymax-ymin);
    ymin = ymin - 0.05 * (ymax-ymin);

    %%-----  INITIAL call  -----
    if k == 0
        %% save default plot bus in the state so we don't have to detect it
        %% each time, since we don't want it to change in the middle of the run
        if plot_bus_default
            cx.cb.default.plot_bus_default = plot_bus_default;
        end
        
        %% initialize lambda-V nose curve plot
        axis([xmin xmax ymin ymax]);
        plot(cxx.lam_hat(1), abs(cxx.V_hat(idx,1)), '-', 'Color', [0.25 0.25 1]);
        title(sprintf('Voltage at Bus %d', idx_e));
        xlabel('\lambda');
        ylabel('Voltage Magnitude');
        hold on;
    %%-----  ITERATION call  -----
    elseif k > 0
        %% plot single step of the lambda-V nose curve
        if plot_level > 1
            axis([xmin xmax ymin ymax]);
            plot([nxx.lam(k); nxx.lam_hat(k+1)], ...
                [abs(nxx.V(idx,k)); abs(nxx.V_hat(idx,k+1))], '-', ...
                'Color', 0.85*[1 0.75 0.75]);
            plot([nxx.lam_hat(k+1); nxx.lam(k+1)], ...
                [abs(nxx.V_hat(idx,k+1)); abs(nxx.V(idx,k+1))], '-', ...
                'Color', 0.85*[0.75 1 0.75]);
            plot(nxx.lam_hat(k+1), abs(nxx.V_hat(idx,k+1)), 'x', ...
                'Color', 0.85*[1 0.75 0.75]);
            plot(nxx.lam(k+1)', abs(nxx.V(idx,k+1))', '-o', ...
                'Color', [0.25 0.25 1]);
            drawnow;
            if plot_level > 2
                pause;
            end
        end
    %%-----  FINAL call  -----
    else    % k < 0
        %% finish final lambda-V nose curve plot
        axis([xmin xmax ymin ymax]);
        plot(nxx.lam', abs(nxx.V(idx,:))',  '-', 'Color', [0.25 0.25 1]);
        hold off;
    end
end
