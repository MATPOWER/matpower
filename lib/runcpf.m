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
%       Alternatively, for compatibility with previous versions of MATPOWER,
%       some of the results can be returned as individual output arguments:
%
%       [baseMVA, bus, gen, branch, success, et] = runcpf(...);
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
%       - if all reference and PV buses are converted, then runcpf raises
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
%   Copyright (c) 1996-2015 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell,
%   Shrirang Abhyankar, Argonne National Laboratory,
%   and Alexander Flueck, IIT
%
%   Modified by Shrirang Abhyankar, Argonne National Laboratory
%   2015.11.25 (Updated to support voltage dependent loads, active and
%               reactive power limits, discrete event handling capabilities)
%
%   $Id$
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

%% define event status codes
cpf_es = cpf_event_status_codes;

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
step             = mpopt.cpf.step;              %% continuation step length
parameterization = mpopt.cpf.parameterization;  %% parameterization
adapt_step       = mpopt.cpf.adapt_step;        %% use adaptive step size?
cb_args          = mpopt.cpf.user_callback_args;
plim             = mpopt.cpf.enforce_p_lims;    %% enforce active limits
qlim             = mpopt.cpf.enforce_q_lims;    %% enforce reactive limits
%% Convert mpopt.cpf.stop_at option to real number. This is to avoid string
%% comparison in event handling
if(ischar(mpopt.cpf.stop_at))
    if(strcmpi(mpopt.cpf.stop_at,'FULL'))
        mpopt.cpf.stop_at = -1;
    else %% NOSE
        mpopt.cpf.stop_at = -2;
    end
end

%% set up callbacks
callback_names = {'cpf_default_callback'};
if ~isempty(mpopt.cpf.user_callback)
    if iscell(mpopt.cpf.user_callback)
        callback_names = {callback_names{:}, mpopt.cpf.user_callback{:}};
    else
        callback_names = {callback_names{:}, mpopt.cpf.user_callback};
    end
end
callbacks = cellfun(@str2func, callback_names, 'UniformOutput', false);

%% run base case power flow
if mpopt.verbose > 2
    mpopt_pf = mpoption(mpopt, 'verbose', max(0, mpopt.verbose-1));
else
    mpopt_pf = mpoption(mpopt, 'verbose', max(0, mpopt.verbose-2));
end
mpopt_pf.pf.enforce_q_lims = mpopt.cpf.enforce_q_lims;
[mpcbase, suc] = runpf(basecasedata, mpopt_pf);
if ~suc
    %% RDZ: this should not be a fatal error
    error('runcpf: Base case power flow did not converge.');
end

%% convert to internal indexing
mpcbase = ext2int(mpcbase);

%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(mpcbase.bus, mpcbase.gen);

nb = size(mpcbase.bus, 1);

%% generator info
onb = find(mpcbase.gen(:, GEN_STATUS) > 0 ...   %% which generators are on?
          & mpcbase.bus(mpcbase.gen(:, GEN_BUS), BUS_TYPE) ~= PQ);  %% ... and not at PQ buses
gbusb = mpcbase.gen(onb, GEN_BUS);             %% what buses are they at?

%% read target case data
mpctarget = loadcase(targetcasedata);

%% add zero columns to branch for flows if needed
if size(mpctarget.branch,2) < QT
  mpctarget.branch = [ mpctarget.branch zeros(size(mpctarget.branch, 1), QT-size(mpctarget.branch,2)) ];
end

%% convert to internal indexing
mpctarget = ext2int(mpctarget);

%% update target bus types in case they've changed due to reactive power limits
mpctarget.bus(:,BUS_TYPE) = mpcbase.bus(:,BUS_TYPE);

%% generator info
ont = find(mpctarget.gen(:, GEN_STATUS) > 0 ...   %% which generators are on?
          & mpctarget.bus(mpctarget.gen(:, GEN_BUS), BUS_TYPE) ~= PQ);  %% ... and not at PQ buses
gbust = mpctarget.gen(ont, GEN_BUS);             %% what buses are they at?

%% RDZ: figure out purpose of this block of code
%% need to ensure the following is same for basecase and targetcase
mpctarget.gen(ont,QG) = mpctarget.gen(onb,QG);
for k = 1:length(ref)
    refgen = find(gbust == ref(k));
    mpctarget.gen(ont(refgen), PG) = mpctarget.gen(onb(refgen), PG);
end

idx_pmax = [];
if plim
    %% Set zero transfer for generators at active power limits
    on = find(mpcbase.gen(:, GEN_STATUS) > 0);
    idx_pmax = find(abs(mpcbase.gen(on,PG) - mpcbase.gen(on,PMAX)) < mpopt.cpf.q_lims_tol);
    if ~isempty(idx_pmax)
        mpctarget.gen(idx_pmax,PG) = mpcbase.gen(idx_pmax,PG);
    end
end

[baseMVAb, busb, genb, branchb] = deal(mpcbase.baseMVA, mpcbase.bus, mpcbase.gen, mpcbase.branch);
[baseMVAt, bust, gent, brancht] = deal(mpctarget.baseMVA, mpctarget.bus, mpctarget.gen, mpctarget.branch);

%%-----  run the continuation power flow  -----
t0 = clock;
if mpopt.verbose > 0
    v = mpver('all');
    fprintf('\nMATPOWER Version %s, %s', v.Version, v.Date);
    fprintf(' -- AC Continuation Power Flow\n');
end
%% initial state
%V0    = ones(size(bus, 1), 1);         %% flat start
V0  = busb(:, VM) .* exp(sqrt(-1) * pi/180 * busb(:, VA));
vcb = ones(size(V0));           %% create mask of voltage-controlled buses
vcb(pq) = 0;                    %% exclude PQ buses
k = find(vcb(gbusb));           %% in-service gens at v-c buses
V0(gbusb(k)) = genb(onb(k), VG) ./ abs(V0(gbusb(k))).* V0(gbusb(k));

%% build admittance matrices
[Ybus, Yf, Yt] = makeYbus(baseMVAb, busb, branchb);

%% function for computing base case V dependent complex bus power injections
%% (generation - load)
Sbusb = @(Vm)makeSbus(baseMVAb, busb, genb, mpopt, Vm);
%% function for computing target case V dependent complex bus power injections
%% (generation - load)
Sbust = @(Vm)makeSbus(baseMVAt, bust, gent, mpopt, Vm);

%% base case power flow solution
lam = 0;
iterations = 0; %% RDZ: modify runpf() so we can get actual from solved power flow
V = V0;
if mpopt.verbose > 2
    fprintf('step %3d : lambda = %6.3f\n', 0, 0);
elseif mpopt.verbose > 1
    fprintf('step %3d : lambda = %6.3f, %2d Newton steps\n', 0, 0, iterations);
end

Vm = abs(V);
lamprv = lam;   %% lam at previous step
Vprv   = V;     %% V at previous step
continuation = 1;
cont_steps = 0;

%% input args for callbacks
Sxfr = Sbust(Vm) - Sbusb(Vm);
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
cb_state = struct();

%% invoke callbacks
for k = 1:length(callbacks)
    cb_state = callbacks{k}(cont_steps, V, lam, V, lam, ...
                            cb_data, cb_state, cb_args);
end

if norm(Sxfr) == 0
    if mpopt.verbose
        fprintf('base case and target case have identical load and generation\n');
    end
    continuation = 0;
end

%% tangent predictor z = [dx;dlam]
z = zeros(2*length(V)+1,1);
z(end,1) = 1.0;
zprv = z;

%% set up events
event_fcn_names = {};
event_postfcn_names = {};
event.names = {};
event.tol = [];
%% Struct for logging information of located events
event.log.stepnums = []; %% Continuation step numbers
event.log.names = {};    %% Names of the events located  
event.log.izero = {};    %% Indices of the events located in the event function
event.log.nevents = [];  %% Number of events located at a step
if qlim
    event_fcn_names = {event_fcn_names{:},'cpf_q_lims_event'};
    event_postfcn_names = {event_postfcn_names{:},'cpf_q_lims_postevent'};
    event.names = {event.names{:},'QGLIM'};
    event.tol = [event.tol,mpopt.cpf.q_lims_tol/baseMVAb];
end
if plim
    event_fcn_names = {event_fcn_names{:},'cpf_p_lims_event'};
    event_postfcn_names = {event_postfcn_names{:},'cpf_p_lims_postevent'};
    event.names = {event.names{:},'PGLIM'};
    event.tol = [event.tol, mpopt.cpf.p_lims_tol/baseMVAb];
end
%% RDZ: order of the above shouldn't matter, but it did in 2015-10-29 version
%%      (linked to tolerances in cpfeventhandler())
%%      not sure if it still matters in 2015-11-25 version
event_fcn_names = {event_fcn_names{:},'cpf_stopat_event'};
event_postfcn_names = {event_postfcn_names{:},'cpf_stopat_postevent'};
if mpopt.cpf.stop_at == -1
    event_name = 'LAM_TRACE_FULL';
elseif mpopt.cpf.stop_at == -2
    event_name = 'LAM_TRACE_NOSE';
else
    event_name = ['LAM_STOPAT',num2str(mpopt.cpf.stop_at)];
end
event.names = {event.names{:},event_name};
event.tol = [event.tol, 1e-5];

event.fcns = cellfun(@str2func,event_fcn_names,'UniformOutput', false);
event.postfcns = cellfun(@str2func,event_postfcn_names,'UniformOutput', false);
event.status = cpf_es.NO_EVENT;
event.qlim_at_prev_step   = 0;

%% prediction for next step
[V0, lam0, z] = cpf_predictor(V, lam, Ybus, Sbusb, Sbust, pv, pq, ...
        step, zprv, Vprv, lamprv, parameterization);
zprv = z;

zprv2 = zprv;
Vprv2 = Vprv;
lamprv2 = lamprv;

%% Call event function
for k = 1:length(event.fcns)
    [event.fprv{k},event.terminate{k}] = event.fcns{k}(cont_steps,V,lam,Vprv,lamprv,z,cb_data);
end

while continuation && event.status ~= cpf_es.TERMINATE
    cont_steps = cont_steps + 1;
    
    %% correction
    [V, success, i, lam] = cpf_corrector(Ybus, Sbusb, V0, ref, pv, pq, ...
                lam0, Sbust, Vprv, lamprv, z, step, parameterization, mpopt_pf);
    if ~success
        continuation = 0;
        if mpopt.verbose
            fprintf('step %3d : lambda = %6.3f, corrector did not converge in %d iterations\n', cont_steps, lam, i);
        end
        break;
    end
    if mpopt.verbose > 2
        fprintf('step %3d : lambda = %6.3f\n', cont_steps, lam);
    elseif mpopt.verbose > 1
        fprintf('step %3d : lambda = %6.3f, %2d corrector Newton steps\n', cont_steps, lam, i);
    end
    
    [V0next, lam0next, z] = cpf_predictor(V, lam, Ybus, Sbusb, Sbust, pv, pq, step, zprv, ...
                              Vprv, lamprv, parameterization);

    %% invoke event handler
    [cont_steps, V, lam, step, cb_data, event, parameterization] = ...
        cpfeventhandler(cont_steps, V, lam, Vprv, lamprv, z, step, ...
            cb_data, event);
    
    if event.status == cpf_es.ZERO_LOC 
        %% Post function may have updated cb_data fields.
        
        ref  = cb_data.ref;
        pv   = cb_data.pv;
        pq   = cb_data.pq;
        
        busb = cb_data.mpc_base.bus;
        genb = cb_data.mpc_base.gen;
        branchb = cb_data.mpc_base.branch;
        
        bust = cb_data.mpc_target.bus;
        gent = cb_data.mpc_target.gen;
        brancht = cb_data.mpc_target.branch;
        
        Sbusb = @(Vm)makeSbus(baseMVAb, busb, genb, mpopt, Vm);
        Sbust = @(Vm)makeSbus(baseMVAt, bust, gent, mpopt, Vm);
        Sxfr  = @(Vm)Sbust(Vm) - Sbusb(Vm);
        
        cb_data.Sbusb = Sbusb;
        cb_data.Sbust = Sbust;
        cb_data.Sxfr  = Sxfr;
        
        [V0next, lam0next, z] = cpf_predictor(V, lam, Ybus, Sbusb, Sbust, pv, pq, step, z, ...
                                              Vprv, lamprv, parameterization);
    end
    
    if ~event.rollback
        %% invoke callbacks
        for k = 1:length(callbacks)
            cb_state = callbacks{k}(cont_steps, V, lam, V0, lam0, ...
                                cb_data, cb_state, cb_args);
        end
        
        V0prv = V0;
        lam0prv = lam0;
        
        V0 = V0next;
        lam0 = lam0next;
        
        %% Update
        Vprv2 = Vprv;
        lamprv2 = lamprv;
        zprv2 = zprv;
        
        Vprv = V;
        lamprv = lam;
        zprv = z;
    else

        %% prediction for next step
        [V0, lam0, z] = cpf_predictor(Vprv, lamprv, Ybus, Sbusb, Sbust, pv, pq, step, zprv2, ...
                                      Vprv2, lamprv2, parameterization);
        
    end
    
    if(~event.rollback && event.status ~= cpf_es.INT_DET)
        if adapt_step && continuation
            %% Adapt stepsize
            cpf_error = norm([angle(V(pq));abs(V([pv;pq]));lam] - [angle(V0prv(pq));abs(V0prv([pv;pq]));lam0prv],inf);
            if cpf_error < mpopt.cpf.error_tol
                %% Increase stepsize
                step = step*mpopt.cpf.error_tol/cpf_error;
                if step > mpopt.cpf.step_max
                    step = mpopt.cpf.step_max;
                end
            else
                %% decrese stepsize
                step = step*mpopt.cpf.error_tol/cpf_error;
                if step < mpopt.cpf.step_min
                    step = mpopt.cpf.step_min;
                end
            end
            
            %% Update prediction
            Vaprv = angle(V);
            Vmprv = abs(V);
            Va0   = angle(V0);
            Vm0   = abs(V0);
            Va0([pv; pq]) = Vaprv([pv; pq]) + step * z([pv; pq]);
            Vm0(pq) = Vmprv(pq) + step * z(nb+pq);
            lam0 = lam + step * z(2*nb+1);
            V0 = Vm0 .* exp(1j * Va0);
        end
    end
end

%% invoke callbacks
if success
    cpf_results = struct();
    for k = 1:length(callbacks)
        [cb_state, cpf_results] = callbacks{k}(cont_steps, V, lam, V0, lam0, ...
                                    cb_data, cb_state, cb_args, cpf_results);                        
    end
else
    cpf_results.iterations = i;
end

%% update bus and gen matrices to reflect the loading and generation
bust(:,PD) = busb(:,PD) + lam*(bust(:,PD) - busb(:,PD));
bust(:,QD) = busb(:,QD) + lam*(bust(:,QD) - busb(:,QD));
gent(:,PG) = genb(:,PG) + lam*(gent(:,PG) - genb(:,PG));

%% update data matrices with solution
[bust, gent, brancht] = pfsoln(baseMVAt, bust, gent, brancht, Ybus, Yf, Yt, V, ref, pv, pq);

mpctarget.et = etime(clock, t0);
mpctarget.success = success;

%%-----  output results  -----
%% convert back to original bus numbering & print results
[mpctarget.bus, mpctarget.gen, mpctarget.branch] = deal(bust, gent, brancht);
if success
    n = cpf_results.iterations + 1;
    cpf_results.V_p = i2e_data(mpctarget, cpf_results.V_p, NaN(nb,n), 'bus', 1);
    cpf_results.V_c = i2e_data(mpctarget, cpf_results.V_c, NaN(nb,n), 'bus', 1);
end
results = int2ext(mpctarget);
results.cpf = cpf_results;
%% Add event log
results.cpf.eventlog = event.log;

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
