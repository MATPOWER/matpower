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
%   If the cpf.enforce_q_lims option is set to true (default is false) then, 
%      -- if any generator reactive power limit is violated during the AC continuation 
%         power flow, the corresponding bus is converted to a PQ bus, with Qg 
%         at the limit. 
%      -- The voltage magnitude at the bus will deviate from the specified value 
%         in order to satisfy the reactive power limit. 
%      -- If the reference bus is converted to PQ, the first remaining PV bus 
%         will be used as the slack bus. This may result in the transfer at this 
%         generator and the generator at the new reference bus being slightly off
%         from the specified values. 
%      -- If all reference and PV buses are converted, then runcpf raises infeasibility 
%         flag and terminates.
%
%   CPF termination modes:
%   Reached nose point 
%   Traced full curve
%   No ref and pv buses remaining
%   Limit induced bifurcation detected
%   
%   Examples:
%       results = runcpf('case9', 'case9target');
%       results = runcpf('case9', 'case9target', ...
%                           mpoption('cpf.adapt_step', 1));
%       results =
%       runcpf('case9','case9target',mpoption('cpf.enforce_q_lims',1));
%       results =
%       runcpf('case9','case9target',mpoption('cpf.stop_at','FULL'));
%
%   See also MPOPTION, RUNPF.

%   MATPOWER
%   Copyright (c) 1996-2015 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell,
%   Shrirang Abhyankar, Argonne National Laboratory,
%   and Alexander Flueck, IIT
%
%   Modified by Shrirang Abhyankar, Argonne National Laboratory
%   2015.10.25 (Updated to support voltage dependent loads, reactive power
%               limits, discrete event handling capabilities)
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

%% Event status flags
[cpf_es] = cpfeventstatus;

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
qlim             = mpopt.cpf.enforce_q_lims;        %% Enforce Q lims
plim             = mpopt.cpf.enforce_p_lims;    %% Enforce P lims

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

%% set up events
event_fcn_names = {'cpf_stopat_event'};
event_postfcn_names = {'cpf_stopat_postevent'};
if qlim
    event_fcn_names = {event_fcn_names{:},'cpf_q_lims_event'};
    event_postfcn_names = {event_postfcn_names{:},'cpf_q_lims_postevent'};
end
if plim
    event_fcn_names = {event_fcn_names{:},'cpf_p_lims_event'};
    event_postfcn_names = {event_postfcn_names{:},'cpf_p_lims_postevent'};
end
event.fcns = cellfun(@str2func,event_fcn_names,'UniformOutput', false);
event.postfcns = cellfun(@str2func,event_postfcn_names,'UniformOutput', false);
event.status = cpf_es.NO_EVENT;

%% run base case power flow
mpopt_pf = mpopt;
mpopt_pf.pf.enforce_q_lims = mpopt.cpf.enforce_q_lims;
mpopt_pf.verbose = mpopt.verbose;
[mpcbase,suc] = runpf(basecasedata,mpopt_pf);
if ~suc
    error('Base case power flow did not converge\n');
end
nb = size(mpcbase.bus, 1);

%% add zero columns to branch for flows if needed
if size(mpcbase.branch,2) < QT
  mpcbase.branch = [ mpcbase.branch zeros(size(mpcbase.branch, 1), QT-size(mpcbase.branch,2)) ];
end

%% convert to internal indexing
mpcbase = ext2int(mpcbase);
[baseMVAb, busb, genb, branchb] = deal(mpcbase.baseMVA, mpcbase.bus, mpcbase.gen, mpcbase.branch);

%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(busb, genb);

%% generator info
%% generator info
onb = find(genb(:, GEN_STATUS) > 0 ...   %% which generators are on?
          & busb(genb(:, GEN_BUS), BUS_TYPE) ~= PQ);  %% ... and not at PQ buses
gbusb = genb(onb, GEN_BUS);             %% what buses are they at?

%% read target case data
mpctarget = loadcase(targetcasedata);

%% add zero columns to branch for flows if needed
if size(mpctarget.branch,2) < QT
  mpctarget.branch = [ mpctarget.branch zeros(size(mpctarget.branch, 1), QT-size(mpctarget.branch,2)) ];
end

%% convert to internal indexing
mpctarget = ext2int(mpctarget);
[baseMVAt, bust, gent, brancht] = deal(mpctarget.baseMVA, mpctarget.bus, mpctarget.gen, mpctarget.branch);
bust(:,BUS_TYPE) = busb(:,BUS_TYPE); %% busb(:,BUS_TYPE) may change if reactive power limits are enforced.

%% generator info
ont = find(gent(:, GEN_STATUS) > 0 ...   %% which generators are on?
          & bust(gent(:, GEN_BUS), BUS_TYPE) ~= PQ);  %% ... and not at PQ buses
gbust = gent(ont, GEN_BUS);             %% what buses are they at?


%% Need to ensure the following is same for basecase and targetcase
gent(ont,QG) = genb(onb,QG);
for k = 1:length(ref)
    refgen = find(gbust == ref(k));
    gent(ont(refgen), PG) = genb(onb(refgen), PG);
end

%%-----  run the power flow  -----
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

%% Run the base case power flow solution
if mpopt.verbose > 2
    mpopt_pf = mpoption(mpopt, 'verbose', max(0, mpopt.verbose-1));
else
    mpopt_pf = mpoption(mpopt, 'verbose', max(0, mpopt.verbose-2));
end
lam = 0;

if mpopt.verbose > 2
    fprintf('step %3d : lambda = %6.3f\n', 0, 0);
elseif mpopt.verbose > 1
    fprintf('step %3d : lambda = %6.3f, %2d Newton steps\n', 0, 0, 0);
end

V = V0;

Sxfr = @(Vm)Sbust(Vm) - Sbusb(Vm);
lamprv = lam;   %% lam at previous step
Vprv   = V;     %% V at previous step
continuation = 1;
cont_steps = 0;

%% input args for callbacks
Sxf = Sxfr(abs(V));
cpf_data = struct( ...
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
    'mpopt', mpopt );
cpf_state = struct();

%% invoke callbacks
for k = 1:length(callbacks)
    cpf_state = callbacks{k}(cont_steps, V, lam, V, lam, ...
                            cpf_data, cpf_state, cb_args);
end

%% event function
for k = 1:length(event.fcns)
    [event.fprv{k},event.terminate{k},cpf_state] = event.fcns{k}(cont_steps,V,lam,cpf_data,cpf_state);
end

if norm(Sxf) == 0
    if mpopt.verbose
        fprintf('base case and target case have identical load and generation\n');
    end
    continuation = 0;
    V0 = V; lam0 = lam;
end

%% tangent predictor z = [dx;dlam]
z = zeros(2*length(V)+1,1);
z(end,1) = 1.0;
zprv = z;
while(continuation && event.status ~= cpf_es.TERMINATE)
    cont_steps = cont_steps + 1;
    Vprv2 = Vprv;
    lamprv2 = lamprv;
    
    %% prediction for next step
    [V0, lam0, z] = cpf_predictor(V, lam, Ybus, Sxfr, Sbust, Sbusb, pv, pq, step, zprv, ...
        Vprv, lamprv, parameterization);

    %% save previous voltage, lambda before updating
    Vprv = V;
    lamprv = lam;
    

    %% correction
    [V, success, i, lam] = cpf_corrector(Ybus, Sbusb, V0, ref, pv, pq, ...
                lam0, Sxfr, Sbust, Vprv, lamprv, z, step, parameterization, mpopt_pf);
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

    %% invoke event handler
    [cont_steps, V, lam, step, cpf_data, cpf_state, event] = cpfeventhandler(cont_steps,V,lam,Vprv,lamprv,step,cpf_data,cpf_state,event);
    
    if event.status == cpf_es.ZERO_LOC 
        %% Post function may have updated cpf_data fields.
        
        ref  = cpf_data.ref;
        pv   = cpf_data.pv;
        pq   = cpf_data.pq;
        
        busb = cpf_data.mpc_base.bus;
        genb = cpf_data.mpc_base.gen;
        branchb = cpf_data.mpc_base.branch;
        
        bust = cpf_data.mpc_target.bus;
        gent = cpf_data.mpc_target.gen;
        brancht = cpf_data.mpc_target.branch;
        
        Sbusb = @(Vm)makeSbus(baseMVAb, busb, genb, mpopt, Vm);
        Sbust = @(Vm)makeSbus(baseMVAt, bust, gent, mpopt, Vm);
        Sxfr  = @(Vm)Sbust(Vm) - Sbusb(Vm);
        
        cpf_data.Sbusb = @(Vm)makeSbus(baseMVAb, busb, genb, mpopt, Vm);
        cpf_data.Sbust = @(Vm)makeSbus(baseMVAt, bust, gent, mpopt, Vm);
        cpf_data.Sxfr  = @(Vm)Sbust(Vm) - Sbusb(Vm);
    end
    
    if ~event.rollback
        %% invoke callbacks
        for k = 1:length(callbacks)
            cpf_state = callbacks{k}(cont_steps, V, lam, V0, lam0, ...
                                cpf_data, cpf_state, cb_args);
        end
        
        %% Update zprv
        zprv = z;
    else
        Vprv = Vprv2;
        lamprv = lamprv2;
    end
    
    if(event.status == cpf_es.INFEASIBLE)
        break;
    end
    
    if(~event.rollback && event.status ~= cpf_es.INT_DET)
        if adapt_step && continuation
            %% Adapt stepsize
            cpf_error = norm([angle(V(pq));abs(V([pv;pq]));lam] - [angle(V0(pq));abs(V0([pv;pq]));lam0],inf);
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
        end
    end
end

%% invoke callbacks
if success
    cpf_results = struct();
    for k = 1:length(callbacks)
        [cpf_state, cpf_results] = callbacks{k}(cont_steps, V, lam, V0, lam0, ...
                                    cpf_data, cpf_state, cb_args, cpf_results);                        
    end
    %% Remove fields added from Qlims and stopat events 
    fieldsrm = {'lam_stopat','lprv','qlim','terminate'};
    cpf_results = rmfield(cpf_results,fieldsrm); 
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

%% Infeasibility 
if(event.status == cpf_es.INFEASIBLE)
    fprintf('CPF ABRUPT TERMINATION!! No ref or PV buses remaining in the system\n');
end

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
