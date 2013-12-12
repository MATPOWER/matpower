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
%   Generator reactive power limits are ignored for the continuation power flow.
%
%   Examples:
%       results = runcpf('case9', 'case9target');
%       results = runcpf('case9', 'case9target', ...
%                           mpoption('cpf.adapt_step', 1));
%
%   See also MPOPTION, RUNPF.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell,
%   Shrirang Abhyankar, Argonne National Laboratory,
%   and Alexander Flueck, IIT
%   Copyright (c) 1996-2013 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

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
step             = mpopt.cpf.step;              %% continuation step length
parameterization = mpopt.cpf.parameterization;  %% parameterization
adapt_step       = mpopt.cpf.adapt_step;        %% use adaptive step size?
cb_args          = mpopt.cpf.user_callback_args;

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

%% read base case data
mpcbase = loadcase(basecasedata);
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
onb = find(genb(:, GEN_STATUS) > 0);      %% which generators are on?
gbusb = genb(onb, GEN_BUS);                %% what buses are they at?

%% read target case data
mpctarget = loadcase(targetcasedata);

%% add zero columns to branch for flows if needed
if size(mpctarget.branch,2) < QT
  mpctarget.branch = [ mpctarget.branch zeros(size(mpctarget.branch, 1), QT-size(mpctarget.branch,2)) ];
end

%% convert to internal indexing
mpctarget = ext2int(mpctarget);
[baseMVAt, bust, gent, brancht] = deal(mpctarget.baseMVA, mpctarget.bus, mpctarget.gen, mpctarget.branch);

%% get bus index lists of each type of bus
%[ref, pv, pq] = bustypes(bust, gent);

%% generator info
ont = find(gent(:, GEN_STATUS) > 0);    %% which generators are on?
gbust = gent(ont, GEN_BUS);             %% what buses are they at?

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

%% compute base case complex bus power injections (generation - load)
Sbusb = makeSbus(baseMVAb, busb, genb);
%% compute target case complex bus power injections (generation - load)
Sbust = makeSbus(baseMVAt, bust, gent);

%% scheduled transfer
Sxfr = Sbust - Sbusb;

%% Run the base case power flow solution
if mpopt.verbose > 2
    mpopt_pf = mpoption(mpopt, 'verbose', max(0, mpopt.verbose-1));
else
    mpopt_pf = mpoption(mpopt, 'verbose', max(0, mpopt.verbose-2));
end
lam = 0;
[V, success, iterations] = newtonpf(Ybus, Sbusb, V0, ref, pv, pq, mpopt_pf);
if mpopt.verbose > 2
    fprintf('step %3d : lambda = %6.3f\n', 0, 0);
elseif mpopt.verbose > 1
    fprintf('step %3d : lambda = %6.3f, %2d Newton steps\n', 0, 0, iterations);
end

lamprv = lam;   %% lam at previous step
Vprv   = V;     %% V at previous step
continuation = 1;
cont_steps = 0;

%% input args for callbacks
cb_data = struct( ...
    'mpc_base', mpcbase, ...
    'mpc_target', mpctarget, ...
    'Sxfr', Sxfr, ...
    'Ybus', Ybus, ...
    'Yf', Yf, ...
    'Yt', Yt, ...
    'ref', ref, ...
    'pv', pv, ...
    'pq', pq, ...
    'mpopt', mpopt );
cb_state = struct();

%% invoke callbacks
for k = 1:length(callbacks)
    cb_state = callbacks{k}(cont_steps, V, lam, V, lam, ...
                            cb_data, cb_state, cb_args);
end

%% tangent predictor z = [dx;dlam]
z = zeros(2*length(V)+1,1);
z(end,1) = 1.0;
while(continuation)
    cont_steps = cont_steps + 1;
    %% prediction for next step
    [V0, lam0, z] = cpf_predictor(V, lam, Ybus, Sxfr, pv, pq, step, z, ...
        Vprv, lamprv, parameterization);

    %% save previous voltage, lambda before updating
    Vprv = V;
    lamprv = lam;

    %% correction
    [V, success, i, lam] = cpf_corrector(Ybus, Sbusb, V0, ref, pv, pq, ...
                lam0, Sxfr, Vprv, lamprv, z, step, parameterization, mpopt_pf);
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

    %% invoke callbacks
    for k = 1:length(callbacks)
        cb_state = callbacks{k}(cont_steps, V, lam, V0, lam0, ...
                                cb_data, cb_state, cb_args);
    end
    
    if ischar(mpopt.cpf.stop_at)
        if strcmp(upper(mpopt.cpf.stop_at), 'FULL')
            if abs(lam) < 1e-8                      %% traced the full continuation curve
                if mpopt.verbose
                    fprintf('\nTraced full continuation curve in %d continuation steps\n',cont_steps);
                end
                continuation = 0;
            elseif lam < lamprv && lam - step < 0   %% next step will overshoot
                step = lam;             %% modify step-size
                parameterization = 1;   %% change to natural parameterization
                adapt_step = 0;         %% disable step-adaptivity
            end
        else    %% == 'NOSE'
            if lam < lamprv                         %% reached the nose point
                if mpopt.verbose
                    fprintf('\nReached steady state loading limit in %d continuation steps\n',cont_steps);
                end
                continuation = 0;
            end
        end
    else
        if lam < lamprv                             %% reached the nose point
            if mpopt.verbose
                fprintf('\nReached steady state loading limit in %d continuation steps\n', cont_steps);
            end
            continuation = 0;
        elseif abs(mpopt.cpf.stop_at - lam) < 1e-8  %% reached desired lambda
            if mpopt.verbose
                fprintf('\nReached desired lambda %3.2f in %d continuation steps\n', ...
                    mpopt.cpf.stop_at, cont_steps);
            end
            continuation = 0;
        elseif lam + step > mpopt.cpf.stop_at   %% will reach desired lambda in next step
            step = mpopt.cpf.stop_at - lam; %% modify step-size
            parameterization = 1;           %% change to natural parameterization
            adapt_step = 0;                 %% disable step-adaptivity
        end
    end
    
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
%% at the noise point
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
        printpf(results, fd, mpopt);
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
