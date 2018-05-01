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
%           can be used to specify the parameterization, output options,
%           termination criteria, and more (see also MPOPTION).
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
%                   user callback functions, where default contains fields:
%               done_msg - string with message describing cause of
%                       continuation termination
%               iterations - number of continuation steps performed
%               lam - (nsteps+1) row vector of lambda values from
%                       correction steps
%               lam_hat - (nsteps+1) row vector of lambda values from
%                       prediction steps
%               max_lam - maximum value of lambda in RESULTS.cpf.lam
%               steps - (nsteps+1) row vector of stepsizes taken at each
%                       continuation step
%               V - (nb x nsteps+1) complex bus voltages from
%                       correction steps
%               V_hat - (nb x nsteps+1) complex bus voltages from
%                       prediction steps
%               events - an array of structs of size nevents with the
%                   following fields:
%                    k - continuation step number at which event was located
%                    name - name of event
%                    idx - index(es) of critical elements in corresponding
%                       event function, e.g. index of generator reaching VAr
%                       limit
%                    msg - descriptive text detailing the event
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
%   If the 'cpf.enforce_q_lims' option is set to true (default is false) then,
%   if any generator reaches its reactive power limits during the AC
%   continuation power flow,
%       - the corresponding bus is converted to a PQ bus, and the problem
%         is modified to eliminate further reactive transfer on this bus
%       - the voltage magnitude at the bus will deviate from the specified
%         setpoint to satisfy the reactive power limit,
%       - if the reference bus is converted to PQ, further real power transfer
%         for the bus is also eliminated, and the first remaining PV bus is
%         selected as the new slack, resulting in the transfers at both
%         reference buses potentially deviating from the specified values
%       - if all reference and PV buses are converted to PQ, RUNCPF terminates
%         with an infeasibility message.
%
%   If the 'cpf.enforce_p_lims' option is set to true (default is fals) then,
%   if any generator reaches its maximum active power limit during the AC
%   continuation power flow,
%       - the problem is modified to eliminate further active transfer by
%         this generator
%       - if the generator was at the reference bus, it is converted to PV
%         and the first remaining PV bus is selected as the new slack.
%
%   If the 'cpf.enforce_v_lims' option is set to true (default is false)
%   then the continuation power flow is set to terminate if any bus voltage
%   magnitude exceeds its minimum or maximum limit.
%
%   If the 'cpf.enforce_flow_lims' option is set to true (default is false)
%   then the continuation power flow is set to terminate if any line MVA
%   flow exceeds its rateA limit.
%
%   Possible CPF termination modes:
%       when cpf.stop_at == 'NOSE'
%           - Reached steady state loading limit
%           - Nose point eliminated by limit induced bifurcation
%       when cpf.stop_at == 'FULL'
%           - Traced full continuation curve
%       when cpf.stop_at == <target_lam_val>
%           - Reached desired lambda
%       when cpf.enforce_p_lims == true
%           - All generators at PMAX
%       when cpf.enforce_q_lims == true
%           - No REF or PV buses remaining
%       when cpf.enforce_v_lims == true
%           - Any bus voltage magnitude limit is reached
%       when cpf.enforce_flow_lims == true
%           - Any branch MVA flow limit is reached
%       other
%           - Base case power flow did not converge
%           - Base and target case have identical load and generation
%           - Corrector did not converge
%           - Could not locate <event_name> event
%           - Too many rollback steps triggered by callbacks
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
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
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
qlim        = mpopt.cpf.enforce_q_lims;    %% enforce reactive limits
plim        = mpopt.cpf.enforce_p_lims;    %% enforce active limits
vlim        = mpopt.cpf.enforce_v_lims;    %% enforce voltage magnitude limits
flim        = mpopt.cpf.enforce_flow_lims; %% enforce branch flow limits

%% register event functions (for event detection)
%% and CPF callback functions (for event handling and other tasks)
cpf_events = [];
cpf_callbacks = [];
%% to handle CPF termination
if ischar(mpopt.cpf.stop_at) && strcmp(mpopt.cpf.stop_at, 'NOSE');
    cpf_events   = cpf_register_event(cpf_events, 'NOSE', 'cpf_nose_event', mpopt.cpf.nose_tol, 1);
    cpf_callbacks = cpf_register_callback(cpf_callbacks, 'cpf_nose_event_cb', 51);
else        %% FULL or target lambda
    cpf_events   = cpf_register_event(cpf_events, 'TARGET_LAM', 'cpf_target_lam_event', mpopt.cpf.target_lam_tol, 1);
    cpf_callbacks = cpf_register_callback(cpf_callbacks, 'cpf_target_lam_event_cb', 50);
end
%% to handle branch flow limits
if flim
    cpf_events = cpf_register_event(cpf_events, 'FLIM', 'cpf_flim_event', mpopt.cpf.flow_lims_tol, 1);
    cpf_callbacks = cpf_register_callback(cpf_callbacks, 'cpf_flim_event_cb', 53);
end
%% to handle voltage limits
if vlim
    cpf_events = cpf_register_event(cpf_events, 'VLIM', 'cpf_vlim_event', mpopt.cpf.v_lims_tol, 1);
    cpf_callbacks = cpf_register_callback(cpf_callbacks, 'cpf_vlim_event_cb', 52);
end
%% to handle reactive power limits
if qlim
    cpf_events = cpf_register_event(cpf_events, 'QLIM', 'cpf_qlim_event', mpopt.cpf.q_lims_tol, 1);
    cpf_callbacks = cpf_register_callback(cpf_callbacks, 'cpf_qlim_event_cb', 41);
end
%% to handle active power limits
if plim
    cpf_events = cpf_register_event(cpf_events, 'PLIM', 'cpf_plim_event', mpopt.cpf.p_lims_tol, 1);
    cpf_callbacks = cpf_register_callback(cpf_callbacks, 'cpf_plim_event_cb', 40);
end
cpf_callbacks = cpf_register_callback(cpf_callbacks, 'cpf_default_callback', 0);
%% user callbacks
if ~isempty(mpopt.cpf.user_callback)
    %% convert to cell array, if necessary
    if iscell(mpopt.cpf.user_callback)
        user_callbacks = mpopt.cpf.user_callback;
    else
        user_callbacks = {mpopt.cpf.user_callback};
    end
    for k = 1:length(user_callbacks)
        %% convert each to struct, if necessary
        if isstruct(user_callbacks{k})
            ucb = user_callbacks{k};
        else
            ucb = struct('fcn', user_callbacks{k});
        end
        %% set default priority, args if not provided
        if ~isfield(ucb, 'priority')
            ucb.priority = [];
        end
        if ~isfield(ucb, 'args')
            ucb.args = [];
        end
        %% register the user callback
        cpf_callbacks = cpf_register_callback(cpf_callbacks, ...
            ucb.fcn, ucb.priority, ucb.args);
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
mpopt_pf = mpoption(mpopt_pf, 'pf.enforce_q_lims', qlim);

%% load base case data
mpcb = loadcase(basecasedata);

%% clip base active generator outputs to PMAX, if necessary
idx_pmax = [];      %% indices of generators clipped at PMAX
if plim
    idx_pmax = find( mpcb.gen(:, GEN_STATUS) > 0 & ...
            mpcb.gen(:, PG) - mpcb.gen(:, PMAX) > -mpopt.cpf.p_lims_tol);
    if mpopt.verbose && ~isempty(idx_pmax)
        fprintf('base case real power output of gen %d reduced from %g to %g MW (PMAX)\n', ...
            [idx_pmax mpcb.gen(idx_pmax, PG) mpcb.gen(idx_pmax, PMAX)]');
    end
    mpcb.gen(idx_pmax, PG) = mpcb.gen(idx_pmax, PMAX);
end

%% run base case power flow
[rb, success] = runpf(mpcb, mpopt_pf);
if success
    done = struct('flag', 0, 'msg', '');
    if qlim
        %% find buses that were converted to PQ or REF by initial power flow
        b2ref = rb.bus(:, BUS_TYPE) == REF & mpcb.bus(:, BUS_TYPE) ~= REF;
        b2pq  = rb.bus(:, BUS_TYPE) == PQ  & mpcb.bus(:, BUS_TYPE) ~= PQ;
    end
    mpcb = rb;      %% update base case with solved power flow
else
    done = struct('flag', 1, 'msg', 'Base case power flow did not converge.');
    results = rb;
    results.cpf = struct();
end

if ~done.flag
    %% read target case data
    mpct = loadcase(targetcasedata);
    if size(mpct.branch,2) < QT     %% add zero columns to branch for flows if needed
      mpct.branch = [ mpct.branch zeros(size(mpct.branch, 1), QT-size(mpct.branch,2)) ];
    end

    %% convert both to internal indexing
    mpcb = ext2int(mpcb, mpopt);
    mpct = ext2int(mpct, mpopt);
    i2e_gen = mpcb.order.gen.i2e;
    if qlim
        b2ref = e2i_data(mpcb, b2ref, 'bus');
        b2pq  = e2i_data(mpcb, b2pq,  'bus');
    end
    nb = size(mpcb.bus, 1);

    %% get bus index lists of each type of bus
    [ref, pv, pq] = bustypes(mpcb.bus, mpcb.gen);

    %% generator info
    %% find generators that are ON and at voltage-controlled buses
    ong = find(mpcb.gen(:, GEN_STATUS) > 0 ...
              & mpcb.bus(mpcb.gen(:, GEN_BUS), BUS_TYPE) ~= PQ);
    gbus = mpcb.gen(ong, GEN_BUS);   %% what buses are they at?

    %% make sure target case is same as base case w.r.t
    %% bus types, GEN_STATUS, Qg and slack Pg
    if qlim
        bb = find(b2ref | b2pq);
        mpct.bus(bb, BUS_TYPE) = mpcb.bus(bb, BUS_TYPE);
    end
    if any(mpcb.bus(:, BUS_TYPE) ~= mpct.bus(:, BUS_TYPE))
        error('runcpf: BUS_TYPE of all buses must be the same in base and target cases');
    end
    if any(mpcb.gen(:, GEN_STATUS) ~= mpct.gen(:, GEN_STATUS))
        error('runcpf: GEN_STATUS of all generators must be the same in base and target cases');
    end
    mpct.gen(ong, QG) = mpcb.gen(ong, QG);
    if qlim
        %% find generators that are ON and at buses that were
        %% converted to PQ by initial power flow
        g2pq  = find(mpcb.gen(:, GEN_STATUS) > 0 & b2pq(mpcb.gen(:, GEN_BUS)));
        mpct.gen(g2pq, QG) = mpcb.gen(g2pq, QG);
    end
    for k = 1:length(ref)
        refgen = find(gbus == ref(k));
        mpct.gen(ong(refgen), PG) = mpcb.gen(ong(refgen), PG);
    end

    %% zero transfers for gens that exceed PMAX limits, if necessary
    if plim
        idx_pmax = find( mpcb.gen(:, GEN_STATUS) > 0 & ...
            mpcb.gen(:, PG)   - mpcb.gen(:, PMAX)   > -mpopt.cpf.p_lims_tol & ...
            mpct.gen(:, PG) - mpct.gen(:, PMAX) > -mpopt.cpf.p_lims_tol);
        if ~isempty(idx_pmax)
            if mpopt.verbose
                fprintf('target case real power output of gen %d reduced from %g to %g MW (PMAX)\n', ...
                    [i2e_gen(idx_pmax) mpct.gen(idx_pmax, PG) mpcb.gen(idx_pmax, PG)]');
            end
            mpct.gen(idx_pmax, PG) = mpcb.gen(idx_pmax, PG);
        end
    end

    %%-----  run the continuation power flow  -----
    t0 = tic;
    if mpopt.verbose
        v = mpver('all');
        fprintf('\nMATPOWER Version %s, %s', v.Version, v.Date);
        fprintf(' -- AC Continuation Power Flow\n');
        if mpopt.verbose > 1
            fprintf('step %3d  :                      lambda = %6.3f, %2d Newton steps\n', 0, 0, mpcb.iterations);
        end
    end

    %% build admittance matrices
    [Ybus, Yf, Yt] = makeYbus(mpcb.baseMVA, mpcb.bus, mpcb.branch);

    %% functions for computing base and target case V-dependent complex bus
    %% power injections: (generation - load)
    Sbusb = @(Vm)makeSbus(mpcb.baseMVA, mpcb.bus, mpcb.gen, mpopt, Vm);
    Sbust = @(Vm)makeSbus(mpct.baseMVA, mpct.bus, mpct.gen, mpopt, Vm);

    %% initialize variables
    cont_steps = 0;
    lam = 0;
    V   = mpcb.bus(:, VM) .* exp(sqrt(-1) * pi/180 * mpcb.bus(:, VA));
    rollback = 0;   %% flag to indicate that a step must be rolled back
    locating = 0;   %% flag to indicate that an event interval was detected,
                    %% but the event has not yet been located
    rb_cnt_ef = 0;  %% counter for rollback steps triggered by event function intervals
    rb_cnt_cb = 0;  %% counter for rollback steps triggered directly by callbacks

    %% initialize tangent predictor: z = [dx;dlam]
    z = [zeros(2*nb, 1); 1];
    direction = 1;
    z = cpf_tangent(V, lam, Ybus, Sbusb, Sbust, pv, pq, ...
                    z, V, lam, parm, direction);

    %% initialize state for current continuation step
    cx = struct(...         %% current state
        'lam_hat', lam, ...         %% predicted lambda
        'V_hat', V, ...             %% predicted V
        'lam', lam, ...             %% corrected lambda
        'V', V, ...                 %% corrected V
        'z', z, ...                 %% normalized tangent predictor
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
        'mpc_base', mpcb, ...       %% MATPOWER case struct - base case
        'mpc_target', mpct, ...     %% MATPOWER case struct - target case
        'Sbusb', Sbusb, ...         %% function for computing base bus inj
        'Sbust', Sbust, ...         %% function for computing target bus inj
        'Ybus', Ybus, ...           %% bus admittance matrix
        'Yf', Yf, ...               %% branch admittance matrix, from end
        'Yt', Yt, ...               %% branch admittance matrix, to end
        'ref', ref, ...             %% vector of ref bus indices
        'pv', pv, ...               %% vector of PV bus indices
        'pq', pq, ...               %% vector of PQ bus indices
        'idx_pmax', idx_pmax, ...   %% vector of gen indices of gens at PMAX
        'mpopt', mpopt );           %% MATPOWER option struct

    %% initialize event function values
    for k = 1:nef
        cx.ef{k} = cpf_events(k).fcn(cb_data, cx);
    end

    %% invoke callbacks - "initialize" context
    for k = 1:ncb
        [nx, cx, done, rollback, evnts, cb_data] = cpf_callbacks(k).fcn( ...
            cont_steps, cx, cx, cx, done, 0, [], cb_data, cpf_callbacks(k).args);
    end

    %% check for case with no transfer
    if norm(Sbust(abs(cx.V)) - Sbusb(abs(cx.V))) == 0
        done.flag = 1;
        done.msg = 'Base case and target case have identical load and generation';
    end

    cont_steps = cont_steps + 1;
    px = cx;    %% initialize state for previous continuation step
    while ~done.flag
        %% initialize next candidate with current state
        nx = cx;

        %% prediction for next step
        [nx.V_hat, nx.lam_hat] = cpf_predictor(cx.V, cx.lam, cx.z, cx.step, cb_data.pv, cb_data.pq);

        %% correction
        [nx.V, success, i, nx.lam] = cpf_corrector(Ybus, cb_data.Sbusb, nx.V_hat, cb_data.ref, cb_data.pv, cb_data.pq, ...
                    nx.lam_hat, cb_data.Sbust, cx.V, cx.lam, cx.z, cx.step, cx.parm, mpopt_pf);
        if ~success     %% corrector failed
            done.flag = 1;
            done.msg = sprintf('Corrector did not converge in %d iterations.', i);
            if mpopt.verbose
                fprintf('step %3d  : stepsize = %-9.3g lambda = %6.3f  corrector did not converge in %d iterations\n', cont_steps, cx.step, nx.lam, i);
            end
            cont_steps = max(cont_steps - 1, 1);    %% go back to last step, but not to 0
            break;
        end

        %% compute new tangent direction, based on a previous state: tx
        if nx.step == 0     %% if this is a re-do step, cx and nx are the same
            tx = px;            %% so use px as the previous state
        else                %% otherwise
            tx = cx;            %% use cx as the previous state
        end
        nx.z = cpf_tangent(nx.V, nx.lam, Ybus, cb_data.Sbusb, cb_data.Sbust, ...
                            cb_data.pv, cb_data.pq, ...
                            cx.z, tx.V, tx.lam, nx.parm, direction);

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
            direction = 1;
        end

        %% invoke callbacks - "iterations" context
        rb = rollback;
        for k = 1:ncb
            [nx, cx, done, rollback, evnts, cb_data] = cpf_callbacks(k).fcn( ...
                cont_steps, nx, cx, px, done, rollback, evnts, cb_data, cpf_callbacks(k).args);
        end
        if ~rb && rollback      %% rollback triggered by callback (vs event function interval)
            rb_cnt_cb = rb_cnt_cb + 1;  %% increment rollback counter for callbacks
            if rb_cnt_cb > 26
                done.flag = 1;
                done.msg = 'Too many rollback steps triggered by callbacks!';
            end
        else
            if ~done.flag && evnts(1).zero
                mpce = cpf_current_mpc(cb_data.mpc_base, cb_data.mpc_target, Ybus, Yf, Yt, cb_data.ref, cb_data.pv, cb_data.pq, nx.V, nx.lam, mpopt);
                J = makeJac(mpce);
                opts.tol = 1e-3;
                opts.it = 2*nb;
                %% retain original direction for buses with negative V-Q
                %% sensitivity, otherwise change direction based on tangent
                %% direction and manifold (eigenvalue)
                if isempty(find(nx.z(nb+pq) > 0, 1))
                    direction = sign(nx.z(end)*min(real(eigs(J,1,'SR',opts))));
                end
            end
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
                done, rollback, evnts, cb_data, cpf_callbacks(k).args, cpf_results);
    end
    cpf_results.events = cx.events;     %% copy eventlog to results

    %% update final case with solution
    mpct = cpf_current_mpc(cb_data.mpc_base, cb_data.mpc_target, Ybus, Yf, Yt, cb_data.ref, cb_data.pv, cb_data.pq, cx.V, cx.lam, mpopt);
    mpct.et = toc(t0);
    mpct.success = success;

    %%-----  output results  -----
    %% convert back to original bus numbering & print results
    n = size(cpf_results.V, 2);
    cpf_results.V_hat = i2e_data(mpct, cpf_results.V_hat, NaN(nb,n), 'bus', 1);
    cpf_results.V     = i2e_data(mpct, cpf_results.V,     NaN(nb,n), 'bus', 1);
    results = int2ext(mpct);
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
