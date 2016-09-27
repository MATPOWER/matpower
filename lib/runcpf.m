function [res, suc] = ...
    runcpf(basecasedata, targetcasedata, mpopt, fname, solvedcase)
%RUNCPF  Runs a full AC continuation power flow
%   [RESULTS, SUCCESS] = RUNCPF(BASECASEDATA, TARGETCASEDATA, ...
%                                   MPOPT, FNAME, SOLVEDCASE)
%
%   Runs a full AC continuation power flow using a normalized tangent
%   predictor and selected parameterization scheme, optionally
%   returning a RESULTS struct and SUCCESS flag. Step size can be
%   fixed or adaptive.
%
%   Inputs (all are optional):
%       BASECASEDATA : either a MATPOWER case struct or a string containing
%           the name of the file with the case data defining the base loading
%           and generation (default is 'case9')
%           (see also CASEFORMAT and LOADCASE)
%       TARGETCASEDATA : either a MATPOWER case struct or a string
%           containing the name of the file with the case data defining the
%           target loading and generation (default is 'case9target')
%       MPOPT : MATPOWER options struct to override default options
%           can be used to specify the solution algorithm, output options
%           termination tolerances, and more (see also MPOPTION).
%       FNAME : name of a file to which the pretty-printed output will
%           be appended
%       SOLVEDCASE : name of file to which the solved case will be saved
%           in MATPOWER case format (M-file will be assumed unless the
%           specified name ends with '.mat')
%
%   Outputs (all are optional):
%       RESULTS : results struct, with the following fields:
%           (all fields from the input MATPOWER case, i.e. bus, branch,
%               gen, etc., but with solved voltages, power flows, etc.)
%           order - info used in external <-> internal data conversion
%           et - elapsed time in seconds
%           success - success flag, 1 = succeeded, 0 = failed
%           cpf - CPF output struct whose content depends on any
%               user callback functions. Default contains fields:
%               V_hat - (nb x nsteps+1) complex bus voltages from
%                       predictor steps
%               lam_hat - (nsteps+1) row vector of lambda values from
%                       predictor steps
%               V - (nb x nsteps+1) complex bus voltages from
%                       corrector steps
%               lam - (nsteps+1) row vector of lambda values from
%                       corrector steps
%               max_lam - maximum value of lambda in lam
%               iterations - number of continuation steps performed
%       SUCCESS : the success flag can additionally be returned as
%           a second output argument
%
%   Calling syntax options:
%       results = runcpf;
%       results = runcpf(basecasedata, targetcasedata);
%       results = runcpf(basecasedata, targetcasedata, mpopt);
%       results = runcpf(basecasedata, targetcasedata, mpopt, fname);
%       results = runcpf(basecasedata, targetcasedata, mpopt, fname, solvedcase)
%       [results, success] = runcpf(...);
%
%   RDZ: Need description of what happens when 'cpf.enforce_p_lims' is true.
%
%   If the 'cpf.enforce_q_lims' option is set to true (default is false) then,
%       - if any generator reactive power limit is violated during the AC
%         continuation power flow, the corresponding bus is converted to a PQ
%         bus, with Qg at the limit.
%       - the voltage magnitude at the bus will deviate from the specified
%         value in order to satisfy the reactive power limit.
%       - if the reference bus is converted to PQ, the first remaining PV bus
%         will be used as the slack bus. This may result in the transfer at
%         this generator and the generator at the new reference bus being
%         slightly off from the specified values.
%       - if all reference and PV buses are converted, then RUNCPF raises
%         infeasibility flag and terminates.
%
%   CPF termination modes:
%       Reached nose point
%       Traced full curve
%       No ref and pv buses remaining
%       Limit induced bifurcation detected
%   
%   Examples:
%       results = runcpf('case9', 'case9target');
%       results = runcpf('case9', 'case9target', ...
%                           mpoption('cpf.adapt_step', 1));
%       results = runcpf('case9', 'case9target', ...
%                           mpoption('cpf.enforce_q_lims', 1));
%       results = runcpf('case9', 'case9target', ...
%                           mpoption('cpf.stop_at', 'FULL'));
%
%   See also MPOPTION, RUNPF.

%   MATPOWER
%   Copyright (c) 1996-2016 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell,
%   Shrirang Abhyankar, Argonne National Laboratory,
%   and Alexander Flueck, IIT
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%-----  initialize  -----
%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;


%% default arguments
if nargin < 5
    solvedcase = '';                %% don't save solved case
    if nargin < 4
        fname = '';                 %% don't print results to a file
        if nargin < 3
            mpopt = mpoption;       %% use default options
            if nargin < 2
                targetcasedata = 'case9target';
                if nargin < 1
                    basecasedata = 'case9'; %% default data file is 'case9.m'
                end
            end
        end
    end
end

%% options
step        = mpopt.cpf.step;              %% continuation step length
parm        = mpopt.cpf.parameterization;  %% parameterization
adapt_step  = mpopt.cpf.adapt_step;        %% use adaptive step size?
cb_args     = mpopt.cpf.user_callback_args;
qlim        = mpopt.cpf.enforce_q_lims;    %% enforce reactive limits
plim        = mpopt.cpf.enforce_p_lims;    %% enforce active limits

%% register event functions (for event detection)
%% and CPF callback functions (for event handling and other tasks)
cpf_events = [];
cpf_callbacks = [];
if ischar(mpopt.cpf.stop_at) && strcmp(mpopt.cpf.stop_at, 'NOSE');
    cpf_events   = cpf_register_event(cpf_events, 'NOSE', 'cpf_nose_event', 1e-5, 1);
    cpf_callbacks = cpf_register_callback(cpf_callbacks, 'cpf_nose_event_cb', 51);
else
    cpf_events   = cpf_register_event(cpf_events, 'TARGET_LAM', 'cpf_target_lam_event', 1e-5, 1);
    cpf_callbacks = cpf_register_callback(cpf_callbacks, 'cpf_target_lam_event_cb', 50);
end
if qlim
    cpf_events = cpf_register_event(cpf_events, 'QLIM', 'cpf_qlim_event', mpopt.cpf.q_lims_tol, 1);
    cpf_callbacks = cpf_register_callback(cpf_callbacks, 'cpf_qlim_event_cb', 41);
end
if plim
    cpf_events = cpf_register_event(cpf_events, 'PLIM', 'cpf_plim_event', mpopt.cpf.p_lims_tol, 1);
    cpf_callbacks = cpf_register_callback(cpf_callbacks, 'cpf_plim_event_cb', 40);
end
cpf_callbacks = cpf_register_callback(cpf_callbacks, 'cpf_default_callback', 0);
if ~isempty(mpopt.cpf.user_callback)
    if iscell(mpopt.cpf.user_callback)
        callback_names = mpopt.cpf.user_callback;
    else
        callback_names = {mpopt.cpf.user_callback};
    end
    for k = 1:length(callback_names)
        cpf_callbacks = cpf_register_callback(cpf_callbacks, callback_names{k}, 10);
    end
end
nef = length(cpf_events);       %% number of event functions registered
ncb = length(cpf_callbacks);    %% number of callback functions registered

%% set power flow options
if mpopt.verbose > 4
    mpopt_pf = mpoption(mpopt, 'verbose', 2);
else
    mpopt_pf = mpoption(mpopt, 'verbose', 0);
end
mpopt_pf = mpoption(mpopt_pf, 'pf.enforce_q_lims', mpopt.cpf.enforce_q_lims);

%% load base case data
mpcbase = loadcase(basecasedata);

%% clip base active generator outputs to PMAX, if necessary
idx_pmax = [];      %% indices of generators clipped at PMAX
if plim
    idx_pmax = find( mpcbase.gen(:, GEN_STATUS) > 0 & ...
            mpcbase.gen(:, PG) - mpcbase.gen(:, PMAX) > -mpopt.cpf.p_lims_tol);
    if mpopt.verbose && ~isempty(idx_pmax)
        fprintf('base case real power output of gen %d reduced from %g to %g MW (PMAX)\n', ...
            [idx_pmax mpcbase.gen(idx_pmax, PG) mpcbase.gen(idx_pmax, PMAX)]');
    end
    mpcbase.gen(idx_pmax, PG) = mpcbase.gen(idx_pmax, PMAX);
end

%% run base case power flow
[mpcbase, suc] = runpf(mpcbase, mpopt_pf);
if suc
    done = struct('flag', 0, 'msg', '');
else
    done = struct('flag', 1, 'msg', 'Base case power flow did not converge.');
    results = mpcbase;
    results.cpf = struct();
end

if ~done.flag
    %% convert to internal indexing
    mpcbase = ext2int(mpcbase);
    nb = size(mpcbase.bus, 1);

    %% get bus index lists of each type of bus
    [ref, pv, pq] = bustypes(mpcbase.bus, mpcbase.gen);

    %% read target case data
    mpctarget = loadcase(targetcasedata);

    %% add zero columns to branch for flows if needed
    if size(mpctarget.branch,2) < QT
      mpctarget.branch = [ mpctarget.branch zeros(size(mpctarget.branch, 1), QT-size(mpctarget.branch,2)) ];
    end

    %% convert to internal indexing
    mpctarget = ext2int(mpctarget);
    i2e_gen = mpctarget.order.gen.i2e;

    %% ensure target case has same bus types as base case
    mpctarget.bus(:, BUS_TYPE) = mpcbase.bus(:, BUS_TYPE);

    %% generator info
    %% find generators that are ON and at voltage-controlled buses
    ong = find(mpcbase.gen(:, GEN_STATUS) > 0 ...
              & mpcbase.bus(mpcbase.gen(:, GEN_BUS), BUS_TYPE) ~= PQ);
    gbus = mpcbase.gen(ong, GEN_BUS);   %% what buses are they at?

    %% make sure target case has same GEN_STATUS
    ont = find(mpctarget.gen(:, GEN_STATUS) > 0 ...
              & mpctarget.bus(mpctarget.gen(:, GEN_BUS), BUS_TYPE) ~= PQ);
    if length(ong) ~= length(ont) || any(ong ~= ont)
        error('runcpf: GEN_STATUS of all generators must be the same in base and target cases');
    end

    %% ensure that Qg and slack Pg for target is same as for base
    %% RDZ: why is this necessary?
    mpctarget.gen(ong, QG) = mpcbase.gen(ong, QG);
    for k = 1:length(ref)
        refgen = find(gbus == ref(k));
        mpctarget.gen(ong(refgen), PG) = mpcbase.gen(ong(refgen), PG);
    end

    %% zero transfers for gens that exceed PMAX limits, if necessary
    if plim
        idx_pmax = find( mpcbase.gen(:, GEN_STATUS) > 0 & ...
            mpcbase.gen(:, PG)   - mpcbase.gen(:, PMAX)   > -mpopt.cpf.p_lims_tol & ...
            mpctarget.gen(:, PG) - mpctarget.gen(:, PMAX) > -mpopt.cpf.p_lims_tol);
        if ~isempty(idx_pmax)
            if mpopt.verbose
                fprintf('target case real power output of gen %d reduced from %g to %g MW (PMAX)\n', ...
                    [i2e_gen(idx_pmax) mpctarget.gen(idx_pmax, PG) mpcbase.gen(idx_pmax, PG)]');
            end
            mpctarget.gen(idx_pmax, PG) = mpcbase.gen(idx_pmax, PG);
        end
    end

    [baseMVAb, busb, genb, branchb] = deal(mpcbase.baseMVA, mpcbase.bus, mpcbase.gen, mpcbase.branch);
    [baseMVAt, bust, gent, brancht] = deal(mpctarget.baseMVA, mpctarget.bus, mpctarget.gen, mpctarget.branch);

    %%-----  run the continuation power flow  -----
    t0 = clock;
    if mpopt.verbose
        v = mpver('all');
        fprintf('\nMATPOWER Version %s, %s', v.Version, v.Date);
        fprintf(' -- AC Continuation Power Flow\n');
    end

    %% build admittance matrices
    [Ybus, Yf, Yt] = makeYbus(baseMVAb, busb, branchb);

    %% functions for computing base and target case V dependent complex bus
    %% power injections: (generation - load)
    Sbusb = @(Vm)makeSbus(baseMVAb, busb, genb, mpopt, Vm);
    Sbust = @(Vm)makeSbus(baseMVAt, bust, gent, mpopt, Vm);

    %% initialize variables
    cont_steps = 0;
    iterations = mpcbase.iterations;
    lam = 0;
    V   = busb(:, VM) .* exp(sqrt(-1) * pi/180 * busb(:, VA));

    %% initialize tangent predictor: z = [dx;dlam]
    z = [zeros(2*nb, 1); 1];
    z = cpf_tangent(V, lam, Ybus, Sbusb, Sbust, pv, pq, z, V, lam, parm);

    %% initialize state for current continuation step
    cx = struct(...         %% current state
        'lam_hat', lam, ...         %% predicted lambda
        'V_hat', V, ...             %% predicted V
        'lam', lam, ...             %% corrected lambda
        'V', V, ...                 %% corrected V
        'z', z, ...                 %% tangent predictor
        'default_step', step, ...   %% default step size
        'default_parm', parm, ...   %% default parameterization
        'this_step', [], ...        %% step size for this step only
        'this_parm', [], ...        %% parameterization for this step only
        'step', step, ...           %% current step size
        'parm', parm, ...           %% current parameterization
        'events', [], ...           %% event log
        'cb', struct(), ...         %% user state, for callbacks (replaces cb_state)
        'ef', {cell(nef, 1)} ...    %% event function values
    );

    %% input args for callbacks
    cb_data = struct( ...
        'mpc_base', mpcbase, ...
        'mpc_target', mpctarget, ...
        'Sbusb', Sbusb, ...
        'Sbust', Sbust, ...
        'Ybus', Ybus, ...
        'Yf', Yf, ...
        'Yt', Yt, ...
        'ref', ref, ...
        'pv', pv, ...
        'pq', pq, ...
        'idx_pmax', idx_pmax, ...
        'mpopt', mpopt );

    %% initialize event function values
    for k = 1:nef
        cx.ef{k} = cpf_events(k).fcn(cb_data, cx);
    end

    if mpopt.verbose > 1
        fprintf('step %3d  :                      lambda = %6.3f, %2d Newton steps\n', 0, 0, iterations);
    end

    %% invoke callbacks - "initialize" context
    for k = 1:ncb
        [nx, cx, done, rollback, evnts, cb_data] = cpf_callbacks(k).fcn( ...
            cont_steps, cx, cx, cx, done, 0, [], cb_data, cb_args);
    end

    %% check for case with no transfer
    if norm(Sbust(abs(cx.V)) - Sbusb(abs(cx.V))) == 0
        done.flag = 1;
        done.msg = 'Base case and target case have identical load and generation';
    end

    rollback = 0;
    locating = 0;
    cont_steps = cont_steps + 1;
    rb_cnt_ef = 0;  %% counter for rollback steps triggered by event function intervals
    rb_cnt_cb = 0;  %% counter for rollback steps triggered directly by callbacks
    px = cx;    %% initialize state for previous continuation step
    while ~done.flag
        %% initialize next candidate with current state
        nx = cx;
    
        %% prediction for next step
        [nx.V_hat, nx.lam_hat] = cpf_predictor(cx.V, cx.lam, cx.z, cx.step, cb_data.pv, cb_data.pq);

        %% correction
        [nx.V, success, i, nx.lam] = cpf_corrector(Ybus, cb_data.Sbusb, nx.V_hat, cb_data.ref, cb_data.pv, cb_data.pq, ...
                    nx.lam_hat, cb_data.Sbust, cx.V, cx.lam, cx.z, cx.step, cx.parm, mpopt_pf);
        if ~success
            done.flag = 1;
            done.msg = sprintf('Corrector did not converge in %d iterations.', i);
            if mpopt.verbose
                fprintf('step %3d  : stepsize = %-9.3g lambda = %6.3f  corrector did not converge in %d iterations\n', cont_steps, cx.step, nx.lam, i);
            end
            cont_steps = cont_steps - 1;
            break;
        end
    
        %% compute new tangent direction, based on a previous state: tx
        if nx.step == 0     %% if this is a re-do step, cx and nx are the same
            tx = px;            %% so use px as the previous state
        else                %% otherwise
            tx = cx;            %% use cx as the previous state
        end
        nx.z = cpf_tangent(nx.V, nx.lam, Ybus, cb_data.Sbusb, cb_data.Sbust, cb_data.pv, cb_data.pq, ...
                                    cx.z, tx.V, tx.lam, nx.parm);

        %% detect events
        for k = 1:nef
            nx.ef{k} = cpf_events(k).fcn(cb_data, nx);      %% update event functions
        end
        [rollback, evnts, nx.ef] = cpf_detect_events(cpf_events, nx.ef, cx.ef, nx.step, mpopt.verbose);

        %% adjust step-size to locate event function zero, if necessary
        if rollback                 %% current step overshot
            %% rollback and initialize next step size based on rollback and previous
            rx = nx;                    %% save state we're rolling back from
            rx.evnts = evnts;           %% and critical event info
            cx.this_step = evnts.step_scale * rx.step;
            cx.this_parm = rx.parm;     %% keep same parameterization as last step
            locating = 1;               %% enter "locating" mode (or stay in it)
            rb_cnt_ef = rb_cnt_ef + 1;  %% increment rollback counter for ef intervals
            if rb_cnt_ef > 26
                done.flag = 1;
                done.msg = sprintf('Could not locate %s event!', evnts.name);
            end
            if mpopt.verbose > 3
                loc_msg = sprintf('OVERSHOOT  : f = [%g, <<%g>>], step <-- %.4g', ...
                            cx.ef{evnts.eidx}(evnts.idx(1)), ...
                            rx.ef{evnts.eidx}(evnts.idx(1)), cx.this_step);
            end
        elseif locating
            if evnts(1).zero        %% found the zero!
                %% step size will be reset to previously used default step size
                locating = 0;           %% exit "locating" mode
                rb_cnt_ef = 0;          %% reset rollback counter for ef intervals
                if mpopt.verbose > 3
                    loc_msg = sprintf('ZERO!      : f = %g, step <-- %.4g', ...
                        nx.ef{rx.evnts.eidx}(rx.evnts.idx(1)), nx.default_step);
                end
            else                    %% prev rollback undershot
                %% initialize next step size based on critical event function
                %% values from prev rollback step and current step
                rx_ef = rx.ef{rx.evnts.eidx}(rx.evnts.idx(1));
                cx_ef = nx.ef{rx.evnts.eidx}(rx.evnts.idx(1));
                step_scale = cx_ef / (cx_ef - rx_ef);
                nx.this_step = step_scale * (rx.step - nx.step);
                rb_cnt_ef = 0;          %% reset rollback counter for ef intervals
                if mpopt.verbose > 3
                    loc_msg = sprintf('UNDERSHOOT : f [<<%g>>, %g], step <-- %.4g', ...
                        cx_ef, rx_ef, nx.this_step);
                end
            end
        else
            loc_msg = '';
        end

        %% invoke callbacks - "iterations" context
        rb = rollback;
        for k = 1:ncb
            [nx, cx, done, rollback, evnts, cb_data] = cpf_callbacks(k).fcn( ...
                cont_steps, nx, cx, px, done, rollback, evnts, cb_data, cb_args);
        end
        if ~rb && rollback      %% rollback triggered by callback (vs event function interval)
            rb_cnt_cb = rb_cnt_cb + 1;  %% increment rollback counter for callbacks
            if rb_cnt_cb > 26
                done.flag = 1;
                done.msg = 'Too many rollback steps triggered by callbacks!';
            end
        else
            rb_cnt_cb = 0;              %% reset rollback counter for callbacks
        end

        %% print iteration information
        if mpopt.verbose > 1
            %% set label for rollback step counter
            if rb_cnt_ef
                sub_step = char('a' + rb_cnt_ef - 1);
            elseif rb_cnt_cb
                sub_step = char('A' + rb_cnt_cb - 1);
            else
                sub_step = ' ';
            end
            if mpopt.verbose > 4
                fprintf('step %3d%s : stepsize = %-9.3g lambda = %6.3f', cont_steps, sub_step, cx.step, nx.lam);
            else
                fprintf('step %3d%s : stepsize = %-9.3g lambda = %6.3f  %2d corrector Newton steps', cont_steps, sub_step, cx.step, nx.lam, i);
            end
            if rollback
                fprintf(' ^ ROLLBACK\n');
            else
                fprintf('\n');
            end
            if mpopt.verbose > 3 && ~isempty(loc_msg)
                fprintf('    LOCATING -- %s\n', loc_msg);
            end
        end

        %% log events
        for k = 1:length(evnts)
            if evnts(k).log
                e = struct( 'k', cont_steps, ...
                            'name', evnts(k).name, ...
                            'idx', evnts(k).idx, ...
                            'msg', evnts(k).msg   );
                if isempty(nx.events)
                    nx.events = e;
                else
                    nx.events(end+1) = e;
                end
            end
            if (mpopt.verbose > 2 && evnts(k).log) || ...
                    (mpopt.verbose > 3 && evnts(k).eidx)
                fprintf('    %s\n', evnts(k).msg);
            end
        end
        
        %% adapt stepsize if requested and not terminating, locating a zero
        %% or re-doing a step after changing the problem data
        if adapt_step && ~done.flag && ~locating && ~evnts(1).zero && nx.step ~= 0
            cpf_error = norm([angle(nx.V(    cb_data.pq)); abs(nx.V(    [cb_data.pv;cb_data.pq])); nx.lam] - ...
                             [angle(nx.V_hat(cb_data.pq)); abs(nx.V_hat([cb_data.pv;cb_data.pq])); nx.lam_hat], inf);

            %% new nominal step size is current size * tol/err, but we reduce
            %% the change from the current size by a damping factor and limit
            %% increases to a factor of 2
            step_scale = min(2, 1 + mpopt.cpf.adapt_step_damping * ...
                            (mpopt.cpf.adapt_step_tol/cpf_error - 1));
            nx.default_step = nx.step * step_scale;

            %% limit step-size
            if nx.default_step > mpopt.cpf.step_max
                nx.default_step = mpopt.cpf.step_max;
            end
            if nx.default_step < mpopt.cpf.step_min
                nx.default_step = mpopt.cpf.step_min;
            end
        end

        %% if this is a normal step
        if ~rollback
            px = cx;    %% save current state before update
            cx = nx;    %% update current state to next candidate
            if ~done.flag
                cont_steps = cont_steps + 1;
            end
        end
    
        %% set step size and parameterization, from one-time or defaults
        if isempty(cx.this_step)
            cx.step = cx.default_step;
        else
            cx.step = cx.this_step;
            cx.this_step = [];      %% disable for next time
        end
        if isempty(cx.this_parm)
            cx.parm = cx.default_parm;
        else
            cx.parm = cx.this_parm;
            cx.this_parm = [];      %% disable for next time
        end
    end

    %% invoke callbacks - "final" context
    cpf_results = struct();     %% initialize results struct
    for k = 1:ncb
        [nx, cx, done, rollback, evnts, cb_data, cpf_results] = ...
            cpf_callbacks(k).fcn(-cont_steps, nx, cx, px, ...
                done, rollback, evnts, cb_data, cb_args, cpf_results);
    end
    cpf_results.events = cx.events;     %% copy eventlog to results

    %% update final case with solution
    mpctarget = cpf_current_mpc(cb_data.mpc_base, cb_data.mpc_target, Ybus, Yf, Yt, cb_data.ref, cb_data.pv, cb_data.pq, cx.V, cx.lam, mpopt);
    mpctarget.et = etime(clock, t0);
    mpctarget.success = success;

    %%-----  output results  -----
    %% convert back to original bus numbering & print results
    n = cpf_results.iterations + 1;
    cpf_results.V_hat = i2e_data(mpctarget, cpf_results.V_hat, NaN(nb,n), 'bus', 1);
    cpf_results.V     = i2e_data(mpctarget, cpf_results.V,     NaN(nb,n), 'bus', 1);
    results = int2ext(mpctarget);
    results.cpf = cpf_results;

    %% zero out result fields of out-of-service gens & branches
    if ~isempty(results.order.gen.status.off)
      results.gen(results.order.gen.status.off, [PG QG]) = 0;
    end
    if ~isempty(results.order.branch.status.off)
      results.branch(results.order.branch.status.off, [PF QF PT QT]) = 0;
    end
end

results.cpf.done_msg = done.msg;
if mpopt.verbose
    fprintf('CPF TERMINATION: %s\n', done.msg);
end

if fname
    [fd, msg] = fopen(fname, 'at');
    if fd == -1
        error(msg);
    else
        if mpopt.out.all == 0
            printpf(results, fd, mpoption(mpopt, 'out.all', -1));
        else
            printpf(results, fd, mpopt);
        end
        fclose(fd);
    end
end
printpf(results, 1, mpopt);

%% save solved case
if solvedcase
    savecase(solvedcase, results);
end

if nargout
    res = results;
    if nargout > 1
        suc = success;
    end
% else  %% don't define res, so it doesn't print anything
end
