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
%               V_p - (nb x nsteps+1) complex bus voltages from
%                       predictor steps
%               lam_p - (nsteps+1) row vector of lambda values from
%                       predictor steps
%               V_c - (nb x nsteps+1) complex bus voltages from
%                       corrector steps
%               lam_c - (nsteps+1) row vector of lambda values from
%                       corrector steps
%               max_lam - maximum value of lambda in lam_c
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


%% register event function callbacks (for event detection)
%% and event handling function callbacks
cpf_events = [];
cpf_handlers = [];
if ischar(mpopt.cpf.stop_at) && strcmp(mpopt.cpf.stop_at, 'NOSE');
    cpf_events   = cpf_register_event(cpf_events, 'NOSE', 'cpf_nose_event', 1e-5, 1);
    cpf_handlers = cpf_register_handler(cpf_handlers, 'NOSE', 'cpf_nose_handler');
else
    cpf_events   = cpf_register_event(cpf_events, 'TARGET_LAM', 'cpf_target_lam_event', 1e-5, 1);
    cpf_handlers = cpf_register_handler(cpf_handlers, 'TARGET_LAM', 'cpf_target_lam_handler');
end
if qlim
    cpf_events = cpf_register_event(cpf_events, 'QLIM', 'cpf_qlim_event', mpopt.cpf.q_lims_tol, 1);
    cpf_handlers = cpf_register_handler(cpf_handlers, 'QLIM', 'cpf_qlim_handler');
end
if plim
    cpf_events = cpf_register_event(cpf_events, 'PLIM', 'cpf_plim_event', mpopt.cpf.p_lims_tol, 1);
    cpf_handlers = cpf_register_handler(cpf_handlers, 'PLIM', 'cpf_plim_handler');
end
cpf_handlers = cpf_register_handler(cpf_handlers, 'DEFAULT', 'cpf_default_handler');
if ~isempty(mpopt.cpf.user_callback)
    if iscell(mpopt.cpf.user_callback)
        callback_names = mpopt.cpf.user_callback;
    else
        callback_names = {mpopt.cpf.user_callback};
    end
    for k = 1:length(callback_names)
        cpf_handlers = cpf_register_handler(cpf_handlers, callback_names{k}, callback_names{k});
    end
end

nef = length(cpf_events);       %% number of event functions registered
neh = length(cpf_handlers);     %% number of event handlers registered

%% register event handling function callbacks

% %% set up callbacks
% callback_names = {'cpf_default_callback'};
% if ~isempty(mpopt.cpf.user_callback)
%     if iscell(mpopt.cpf.user_callback)
%         callback_names = {callback_names{:}, mpopt.cpf.user_callback{:}};
%     else
%         callback_names = {callback_names{:}, mpopt.cpf.user_callback};
%     end
% end
% callbacks = cellfun(@str2func, callback_names, 'UniformOutput', false);

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
if ~suc
    %% RDZ: this should not be a fatal error, simply return suc = 0
    %%      possibly by setting continuation to 0, rather than 1, later
    error('runcpf: Base case power flow did not converge.');
end

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
Sxfr  = @(Vm)(Sbust(Vm) - Sbusb(Vm));

%% initialize variables
continuation = 1;
cont_steps = 0;
iterations = mpcbase.iterations;
lam = 0;
V   = busb(:, VM) .* exp(sqrt(-1) * pi/180 * busb(:, VA));

%% initialize tangent predictor: z = [dx;dlam]
z = [zeros(2*nb, 1); 1];
z = cpf_tangent(V, lam, Ybus, Sbusb, Sbust, pv, pq, ...
                            z, V, lam, parm);

%% initialize values for current continuation step
cc = struct(...         %% current values
    'lam0', lam, ...            %% predicted lambda
    'V0', V, ...                %% predicted V
    'lam', lam, ...             %% corrected lambda
    'V', V, ...                 %% corrected V
    'z', z, ...                 %% tangent predictor
    'default_step', step, ...   %% default step size
    'default_parm', parm, ...   %% default parameterization
    'this_step', [], ...        %% step size for this step only
    'this_parm', [], ...        %% parameterization for this step only
    'step', step, ...           %% current step size
    'parm', parm, ...           %% current parameterization
    'ef', {cell(nef, 1)} ...    %% event function values
);

%% input args for callbacks
cb_data = struct( ...
    'mpc_base', mpcbase, ...
    'mpc_target', mpctarget, ...
    'Sbusb', Sbusb, ...
    'Sbust', Sbust, ...
    'Sxfr', Sxfr, ...
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
    cc.ef{k} = cpf_events(k).fcn(cb_data, cc);
end

if mpopt.verbose > 1
    fprintf('step %3d  :                      lambda = %6.3f, %2d Newton steps\n', 0, 0, iterations);
end

%% initialize callback state
cb_state = struct();

%% invoke event handlers - "initialize" context
for k = 1:neh
    [cb_state, nn, cc, cb_data, terminate] = cpf_handlers(k).fcn(cont_steps, ...
        cc, cc, cc, 0, [], 0, cb_data, cb_state, cb_args);
end

% %% invoke callbacks - "initial" context
% for k = 1:length(callbacks)
%     cb_state = callbacks{k}(cont_steps, cc.step, cc.V, cc.lam, cc.V0, cc.lam0, ...
%                             cb_data, cb_state, cb_args);
% end

%% check for case with no transfer
if norm(Sxfr(abs(cc.V))) == 0
    if mpopt.verbose
        fprintf('base case and target case have identical load and generation\n');
    end
    continuation = 0;
end

rollback = 0;
locating = 0;
prev_rollback = rollback;
cont_steps = cont_steps + 1;
sub_step = ' ';
pp = cc;    %% initialize values for previous continuation step
while continuation
    %% initialize next candidate with current values
    nn = cc;
    
    %% prediction for next step
    [nn.V0, nn.lam0] = cpf_predictor(cc.V, cc.lam, cc.z, cc.step, cb_data.pv, cb_data.pq);

    %% correction
    [nn.V, success, i, nn.lam] = cpf_corrector(Ybus, cb_data.Sbusb, nn.V0, cb_data.ref, cb_data.pv, cb_data.pq, ...
                nn.lam0, cb_data.Sbust, cc.V, cc.lam, cc.z, cc.step, cc.parm, mpopt_pf);
    if ~success
        continuation = 0;
        break;
    end
    
    %% compute new tangent direction
    if nn.step == 0
        pV = pp.V;
        plam = pp.lam;
    else
        pV = cc.V;
        plam = cc.lam;
    end
    nn.z = cpf_tangent(nn.V, nn.lam, Ybus, cb_data.Sbusb, cb_data.Sbust, cb_data.pv, cb_data.pq, ...
                                cc.z, pV, plam, nn.parm);

    %% update event functions
    for k = 1:nef
        nn.ef{k} = cpf_events(k).fcn(cb_data, nn);
    end
    
    %% detect events
    prev_rollback = rollback;   %% save rollback flag value from prev step
%mpopt.verbose = 3;
    [rollback, critical, nn.ef] = cpf_detect_events(cpf_events, nn.ef, cc.ef, nn.step, mpopt.verbose);

    %% adjust step-size to locate event function zero, if necessary
    if rollback                                     %% current step overshot
        %% rollback and initialize next step size based on rollback and previous
        rb = nn;                    %% save rolled back values
        rb.critical = critical;     %% and critical event info
        cc.this_step = critical.step_scale * rb.step;
        cc.this_parm = rb.parm;     %% keep same parameterization as last step
        if locating
            sub_step = char(sub_step + 1);
            if sub_step > 'z'
                if mpopt.verbose
                    fprintf('CPF Termination : Could not locate %s event!\n', critical.name);
                end
                continuation = 0;
            end
        else
            locating = 1;           %% enter "locating" mode
            sub_step = 'a';
        end
        if mpopt.verbose > 3
            fprintf('   -- OVERSHOOT  : f = [%g, <<%g>>], step = %g\n', ...
                        cc.ef{critical.k}(critical.idx(1)), ...
                        rb.ef{critical.k}(critical.idx(1)), cc.this_step);
        end
    elseif locating
        if strcmp(critical(1).status, 'ZERO')       %% found the zero!
            %% reset to the previously used default step size
%            nn.step = nn.default_step;
            locating = 0;           %% exit "locating" mode
            sub_step = ' ';
            if mpopt.verbose > 3
                fprintf('   -- ZERO!      : f = %g, step = %g\n', ...
                    nn.ef{rb.critical.k}(rb.critical.idx(1)), nn.default_step);
            end
        else                                        %% prev rollback undershot
            %% initialize next step size based on critical event function
            %% values from prev rollback step and current step
            rbef = rb.ef{rb.critical.k}(rb.critical.idx(1));
            ccef = nn.ef{rb.critical.k}(rb.critical.idx(1));
            step_scale = ccef / (ccef - rbef);
            nn.this_step = step_scale * (rb.step - nn.step);
            sub_step = ' ';
            if mpopt.verbose > 3
                fprintf('   -- UNDERSHOOT : f [<<%g>>, %g], step = %g\n', ccef, rbef, nn.this_step);
            end
        end
%     else
%         if mpopt.verbose > 3
%             fprintf('   -- NORMAL\n');
%         end
    end

    %% invoke event handlers - "iterations" context
    terminate = ~continuation;
    for k = 1:neh
        [cb_state, nn, cc, cb_data, terminate] = cpf_handlers(k).fcn(cont_steps, ...
            nn, cc, pp, rollback, critical, terminate, cb_data, cb_state, cb_args);
        if terminate
            continuation = 0;
        end
    end

    if mpopt.verbose > 4
        fprintf('step %3d%s : stepsize = %-9.3g lambda = %6.3f\n', cont_steps, sub_step, cc.step, nn.lam);
    elseif mpopt.verbose > 1
        fprintf('step %3d%s : stepsize = %-9.3g lambda = %6.3f  %2d corrector Newton steps\n', cont_steps, sub_step, cc.step, nn.lam, i);
    end

% if ~rollback
%     %% invoke callbacks - "iterations" context
%     for k = 1:length(callbacks)
%         cb_state = callbacks{k}(cont_steps, nn.step, nn.V, nn.lam, nn.V0, nn.lam0, ...
%                             cb_data, cb_state, cb_args);
%     end
% end
% 
%     if ischar(mpopt.cpf.stop_at)
%         if strcmp(upper(mpopt.cpf.stop_at), 'FULL')
%             if abs(nn.lam) < 1e-8                   %% traced the full continuation curve
% %             if nn.lam < 1e-8                        %% traced the full continuation curve
%                 if mpopt.verbose
%                     fprintf('\nTraced full continuation curve in %d continuation steps\n',cont_steps);
%                 end
%                 continuation = 0;
%             elseif nn.lam < cc.lam && nn.lam - cc.step < 0  %% next step will overshoot
%                 nn.this_step = nn.lam;  %% modify step-size
%                 nn.this_parm = 1;       %% change to natural parameterization
%                 adapt_step = 0;         %% disable step-adaptivity
%             end
%         else    %% == 'NOSE'
%             if nn.lam < cc.lam                      %% reached the nose point
%                 if mpopt.verbose
%                     fprintf('\nReached steady state loading limit in %d continuation steps\n',cont_steps);
%                 end
%                 continuation = 0;
%             end
%         end
%     else
%         if nn.lam < cc.lam                          %% reached the nose point
%             if mpopt.verbose
%                 fprintf('\nReached steady state loading limit in %d continuation steps\n', cont_steps);
%             end
%             continuation = 0;
%         elseif abs(mpopt.cpf.stop_at - nn.lam) < 1e-8   %% reached desired lambda
%             if mpopt.verbose
%                 fprintf('\nReached desired lambda %3.2f in %d continuation steps\n', ...
%                     mpopt.cpf.stop_at, cont_steps);
%             end
%             continuation = 0;
%         elseif nn.lam + cc.step > mpopt.cpf.stop_at %% will reach desired lambda in next step
%             nn.this_step = mpopt.cpf.stop_at - nn.lam;  %% modify step-size
%             nn.this_parm = 1;               %% change to natural parameterization
%             adapt_step = 0;                 %% disable step-adaptivity
%         end
%     end
    
    if adapt_step && continuation && ~locating && ~strcmp(critical(1).status, 'ZERO') && nn.step ~= 0
        %% adapt stepsize
        cpf_error = norm([angle(nn.V(cb_data.pq));  abs(nn.V([cb_data.pv;cb_data.pq]));  nn.lam] - ...
                         [angle(nn.V0(cb_data.pq)); abs(nn.V0([cb_data.pv;cb_data.pq])); nn.lam0], inf);

        %% new nominal step size is current size * tol/err, but we reduce
        %% the change from the current size by a damping factor and limit
        %% increases to a factor of 2
        %% RDZ: reduce this to 2 for release (and update tests)
        ff = 10000;
        step_scale = min(ff, 1 + mpopt.cpf.adapt_step_damping * ...
                        (mpopt.cpf.error_tol/cpf_error - 1));
        nn.default_step = nn.step * step_scale;

        %% limit step-size
        if nn.default_step > mpopt.cpf.step_max
            nn.default_step = mpopt.cpf.step_max;
        end
        if nn.default_step < mpopt.cpf.step_min
            nn.default_step = mpopt.cpf.step_min;
        end
%fprintf('---- ADAPT ');
    end

    if ~rollback
        pp = cc;    %% save current values before update
        cc = nn;    %% update current point to next candidate
        if continuation
            cont_steps = cont_steps + 1;
        end
    end
    
    %% set current step size and parameterization, from one-time or defaults
    if isempty(cc.this_step)
        cc.step = cc.default_step;
%fprintf('---- DEFAULT : %g\n', cc.step);
    else
%         if cc.this_step == 0    %% this is a "repeat" step (e.g. after bus type changes)
%             rollback = 0;       %% don't treat as a rollback step
%             locating = 0;       %% and exit "locating" mode, too
%         end
        cc.step = cc.this_step;
        cc.this_step = [];      %% disable for next time
%fprintf('---- ONETIME : %g\n', cc.step);
    end
    if isempty(cc.this_parm)
        cc.parm = cc.default_parm;
    else
        cc.parm = cc.this_parm;
        cc.this_parm = [];      %% disable for next time
    end
end

%% invoke callbacks - "final" context
if success
    cpf_results = struct();

    %% invoke event handlers - "finalize" context
    for k = 1:neh
        [cb_state, nn, cc, cb_data, terminate, cpf_results] = cpf_handlers(k).fcn(-cont_steps, ...
            nn, cc, pp, rollback, critical, 0, cb_data, cb_state, cb_args, cpf_results);
    end

%     for k = 1:length(callbacks)
%         [cb_state, cpf_results] = callbacks{k}(cont_steps, cc.step, cc.V, cc.lam, cc.V0, cc.lam0, ...
%                                     cb_data, cb_state, cb_args, cpf_results);
%     end
else
    if mpopt.verbose
        fprintf('step %3d%s : stepsize = %-9.3g lambda = %6.3f  corrector did not converge in %d iterations\n', cont_steps, sub_step, cc.step, nn.lam, i);
    end

    cpf_results.iterations = i;
end

% %% update bus and gen matrices to reflect the loading and generation
% bust(:,PD) = busb(:,PD) + cc.lam*(bust(:,PD) - busb(:,PD));
% bust(:,QD) = busb(:,QD) + cc.lam*(bust(:,QD) - busb(:,QD));
% gent(:,PG) = genb(:,PG) + cc.lam*(gent(:,PG) - genb(:,PG));
% 
% %% update data matrices with solution
% [bust, gent, brancht] = pfsoln(baseMVAt, bust, gent, brancht, Ybus, Yf, Yt, cc.V, cb_data.ref, cb_data.pv, cb_data.pq, mpopt);

%% update final case with solution
mpctarget = cpf_current_mpc(cb_data.mpc_base, cb_data.mpc_target, Ybus, Yf, Yt, cb_data.ref, cb_data.pv, cb_data.pq, cc.V, cc.lam, mpopt);
mpctarget.et = etime(clock, t0);
mpctarget.success = success;

%%-----  output results  -----
%% convert back to original bus numbering & print results
% [mpctarget.bus, mpctarget.gen, mpctarget.branch] = deal(bust, gent, brancht);
if success
    n = cpf_results.iterations + 1;
    cpf_results.V_p = i2e_data(mpctarget, cpf_results.V_p, NaN(nb,n), 'bus', 1);
    cpf_results.V_c = i2e_data(mpctarget, cpf_results.V_c, NaN(nb,n), 'bus', 1);
end
results = int2ext(mpctarget);
results.cpf = cpf_results;

%% zero out result fields of out-of-service gens & branches
if ~isempty(results.order.gen.status.off)
  results.gen(results.order.gen.status.off, [PG QG]) = 0;
end
if ~isempty(results.order.branch.status.off)
  results.branch(results.order.branch.status.off, [PF QF PT QT]) = 0;
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
