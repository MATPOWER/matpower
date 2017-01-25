function [nx, cx, done, rollback, evnts, cb_data, results] = ...
    cpf_default_callback(k, nx, cx, px, done, rollback, evnts, ...
                            cb_data, cb_args, results)
%CPF_DEFAULT_CALLBACK   Default callback function for CPF
%   [NX, CX, DONE, ROLLBACK, EVNTS, CB_DATA, RESULTS] = 
%       CPF_DEFAULT_CALLBACK(K, NX, CX, PX, DONE, ROLLBACK, EVNTS, ...
%                               CB_DATA, CB_ARGS, RESULTS)
%
%   Default callback function used by RUNCPF that collects the resulst and
%   optionally, plots the nose curve. Inputs and outputs are defined below,
%   with the RESULTS argument being optional, used only for the final call
%   when K is negative.
%
%   Inputs:
%       K - continuation step iteration count
%       NX - next state (corresponding to proposed next step), struct with
%            the following fields:
%           lam_hat - value of LAMBDA from predictor
%           V_hat - vector of complex bus voltages from predictor
%           lam - value of LAMBDA from corrector
%           V - vector of complex bus voltages from corrector
%           z - normalized tangent predictor
%           default_step - default step size
%           default_parm - default parameterization
%           this_step - step size for this step only
%           this_parm - paramterization for this step only
%           step - current step size
%           parm - current parameterization
%           events - struct array, event log
%           cb - user state, for callbacks (replaces CB_STATE), the user may
%               add fields containing any information the callback function
%               would like to pass from one invokation to the next, taking
%               care not to step on fields being used by other callbacks,
%               such as the 'default' field used by this default callback
%           ef - cell array of event function values
%       CX - current state (corresponding to most recent successful step)
%            (same structure as NX)
%       PX - previous state (corresponding to last successful step prior to CX)
%       DONE - struct, with flag to indicate CPF termination and reason,
%           with fields:
%           flag - termination flag, 1 => terminate, 0 => continue
%           msg - string containing reason for termination
%       ROLLBACK - scalar flag to indicate that the current step should be
%           rolled back and retried with a different step size, etc.
%       EVNTS - struct array listing any events detected for this step,
%           see CPF_DETECT_EVENTS for details
%       CB_DATA - struct containing potentially useful "static" data,
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
%           idx_pmax - vector of generator indices for generators fixed
%               at their PMAX limits
%           mpopt - MATPOWER options struct
%       CB_ARGS - arbitrary data structure containing callback arguments
%       RESULTS - initial value of output struct to be assigned to
%           CPF field of results struct returned by RUNCPF
%
%   Outputs:
%       (all are updated versions of the corresponding input arguments)
%       NX - user state ('cb' field ) should be updated here if ROLLBACK
%           is false
%       CX - may contain updated 'this_step' or 'this_parm' values to be used
%           if ROLLBACK is true
%       DONE - callback may have requested termination and set the msg field
%       ROLLBACK - callback can request a rollback step, even if it was not
%           indicated by an event function
%       EVNTS - msg field for a given event may be updated
%       CB_DATA - this data should only be modified if the underlying problem
%           has been changed (e.g. generator limit reached) and should always
%           be followed by a step of zero length, i.e. set NX.this_step to 0
%           It is the job of any callback modifying CB_DATA to ensure that
%           all data in CB_DATA is kept consistent.
%       RESULTS - updated version of RESULTS input arg
%
%   This function is called in three different contexts, distinguished by
%   the value of K, as follows:
%   (1) initial - called with K = 0, without RESULTS input/output args,
%           after base power flow, before 1st CPF step.
%   (2) iterations - called with K > 0, without RESULTS input/output args,
%           at each iteration, after predictor-corrector step
%   (3) final - called with K < 0, with RESULTS input/output args, after
%           exiting predictor-corrector loop, inputs identical to last
%           iteration call, except K which is negated
%
%   User Defined CPF Callback Functions:
%       The user can define their own callback functions which take
%       the same form and are called in the same contexts as
%       CPF_DEFAULT_CALLBACK. These are specified via the MATPOWER
%       option 'cpf.user_callback'. This option can be a string containing
%       the name of the callback function, or a struct with the following
%       fields, where all but the first are optional:
%           'fcn'       - string with name of callback function
%           'priority'  - numerical value specifying callback priority
%                (default = 20, see CPF_REGISTER_CALLBACK for details)
%           'args'      - arbitrary value (any type) passed to the callback
%                         as CB_ARGS each time it is invoked
%       Multiple user callbacks can be registered by assigning a cell array
%       of such strings and/or structs to the 'cpf.user_callback' option.
%
%   See also RUNCPF, CPF_REGISTER_CALLBACK.

%   MATPOWER
%   Copyright (c) 2013-2016, Power Systems Engineering Research Center (PSERC)
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
        if any(idx_e > length(cb_data.mpc_target.order.bus.e2i))
            kk = find(idx_e > length(cb_data.mpc_target.order.bus.e2i));
            error('cpf_default_callback: %d is not a valid bus number for MPOPT.cpf.plot.bus', idx_e(kk(1)));
        end
        idx = full(cb_data.mpc_target.order.bus.e2i(idx_e));
        if any(idx == 0)
            kk = find(idx == 0);
            error('cpf_default_callback: %d is not a valid bus number for MPOPT.cpf.plot.bus', idx_e(kk(1)));
        end
    end
    nplots = length(idx_e);

    %% set bounds for plot axes
    xmin = 0;
    xmax = max([max(nxx.lam_hat); max(nxx.lam)]);
    ymin = min([min(min(abs(nxx.V_hat(idx, :)))); min(min(abs(nxx.V(idx, :))))]);
    ymax = max([max(max(abs(nxx.V_hat(idx, :)))); max(max(abs(nxx.V(idx, :))))]);
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
        if nplots > 1
            title('Voltage at Multiple Buses');
        else
            title(sprintf('Voltage at Bus %d', idx_e));
        end
        xlabel('\lambda');
        ylabel('Voltage Magnitude');
        hold on;
    %%-----  ITERATION call  -----
    elseif k > 0
        %% plot single step of the lambda-V nose curve
        if plot_level > 1
            axis([xmin xmax ymin ymax]);
            for kk = 1:nplots
                %% prediction line
                plot([nxx.lam(k); nxx.lam_hat(k+1)], ...
                    [abs(nxx.V(idx(kk),k)); abs(nxx.V_hat(idx(kk),k+1))], '-', ...
                    'Color', 0.85*[1 0.75 0.75]);
                %% correction line
                plot([nxx.lam_hat(k+1); nxx.lam(k+1)], ...
                    [abs(nxx.V_hat(idx(kk),k+1)); abs(nxx.V(idx(kk),k+1))], '-', ...
                    'Color', 0.85*[0.75 1 0.75]);
                %% prediciton point
                plot(nxx.lam_hat(k+1), abs(nxx.V_hat(idx(kk),k+1)), 'x', ...
                    'Color', 0.85*[1 0.75 0.75]);
                %% correction point
                plot(nxx.lam(k+1)', abs(nxx.V(idx(kk),k+1))', '-o', ...
                    'Color', [0.25 0.25 1]);
                drawnow;
            end
            if plot_level > 2
                pause;
            end
        end
    %%-----  FINAL call  -----
    else    % k < 0
        %% finish final lambda-V nose curve plot
        axis([xmin xmax ymin ymax]);
        %% curve of corrected points
        if isprop(gca, 'ColorOrderIndex')
            set(gca, 'ColorOrderIndex', 1); %% start over with color 1
        end
        hp = plot(nxx.lam', abs(nxx.V(idx,:))',  '-');
        if nplots > 1
            leg = cell(nplots, 1);
            for kk = 1:nplots
                leg{kk} = sprintf('Bus %d', idx_e(kk));
            end
            legend(hp, leg);
        end
        hold off;
    end
end
