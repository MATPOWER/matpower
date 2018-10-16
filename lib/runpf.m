function [MVAbase, bus, gen, branch, success, et] = ...
                runpf(casedata, mpopt, fname, solvedcase)
%RUNPF  Runs a power flow.
%   [RESULTS, SUCCESS] = RUNPF(CASEDATA, MPOPT, FNAME, SOLVEDCASE)
%
%   Runs a power flow (full AC Newton's method by default), optionally
%   returning a RESULTS struct and SUCCESS flag.
%
%   Inputs (all are optional):
%       CASEDATA : either a MATPOWER case struct or a string containing
%           the name of the file with the case data (default is 'case9')
%           (see also CASEFORMAT and LOADCASE)
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
%       SUCCESS : the success flag can additionally be returned as
%           a second output argument
%
%   Calling syntax options:
%       results = runpf;
%       results = runpf(casedata);
%       results = runpf(casedata, mpopt);
%       results = runpf(casedata, mpopt, fname);
%       results = runpf(casedata, mpopt, fname, solvedcase);
%       [results, success] = runpf(...);
%
%       Alternatively, for compatibility with previous versions of MATPOWER,
%       some of the results can be returned as individual output arguments:
%
%       [baseMVA, bus, gen, branch, success, et] = runpf(...);
%
%   If the pf.enforce_q_lims option is set to true (default is false) then, if
%   any generator reactive power limit is violated after running the AC power
%   flow, the corresponding bus is converted to a PQ bus, with Qg at the
%   limit, and the case is re-run. The voltage magnitude at the bus will
%   deviate from the specified value in order to satisfy the reactive power
%   limit. If the reference bus is converted to PQ, the first remaining PV
%   bus will be used as the slack bus for the next iteration. This may
%   result in the real power output at this generator being slightly off
%   from the specified values.
%
%   Examples:
%       results = runpf('case30');
%       results = runpf('case30', mpoption('pf.enforce_q_lims', 1));
%
%   See also RUNDCPF.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   Enforcing of generator Q limits inspired by contributions
%   from Mu Lin, Lincoln University, New Zealand (1/14/05).
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
if nargin < 4
    solvedcase = '';                %% don't save solved case
    if nargin < 3
        fname = '';                 %% don't print results to a file
        if nargin < 2
            mpopt = mpoption;       %% use default options
            if nargin < 1
                casedata = 'case9'; %% default data file is 'case9.m'
            end
        end
    end
end

%% options
qlim = mpopt.pf.enforce_q_lims;         %% enforce Q limits on gens?
dc = strcmp(upper(mpopt.model), 'DC');  %% use DC formulation?

%% read data
mpc = loadcase(casedata);

%% add zero columns to branch for flows if needed
if size(mpc.branch,2) < QT
  mpc.branch = [ mpc.branch zeros(size(mpc.branch, 1), QT-size(mpc.branch,2)) ];
end

%% convert to internal indexing
mpc = ext2int(mpc, mpopt);
[baseMVA, bus, gen, branch] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);

if ~isempty(mpc.bus)
    %% get bus index lists of each type of bus
    [ref, pv, pq] = bustypes(bus, gen);

    %% generator info
    on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on?
    gbus = gen(on, GEN_BUS);                %% what buses are they at?

    %%-----  run the power flow  -----
    t0 = tic;
    its = 0;            %% total iterations
    if mpopt.verbose > 0
        v = mpver('all');
        fprintf('\nMATPOWER Version %s, %s', v.Version, v.Date);
    end
    if dc                               %% DC formulation
        if mpopt.verbose > 0
          fprintf(' -- DC Power Flow\n');
        end
        %% initial state
        Va0 = bus(:, VA) * (pi/180);
    
        %% build B matrices and phase shift injections
        [B, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);
    
        %% compute complex bus power injections (generation - load)
        %% adjusted for phase shifters and real shunts
        Pbus = real(makeSbus(baseMVA, bus, gen)) - Pbusinj - bus(:, GS) / baseMVA;
    
        %% "run" the power flow
        [Va, success] = dcpf(B, Pbus, Va0, ref, pv, pq);
        its = 1;
    
        %% update data matrices with solution
        branch(:, [QF, QT]) = zeros(size(branch, 1), 2);
        branch(:, PF) = (Bf * Va + Pfinj) * baseMVA;
        branch(:, PT) = -branch(:, PF);
        bus(:, VM) = ones(size(bus, 1), 1);
        bus(:, VA) = Va * (180/pi);
        %% update Pg for slack generator (1st gen at ref bus)
        %% (note: other gens at ref bus are accounted for in Pbus)
        %%      Pg = Pinj + Pload + Gs
        %%      newPg = oldPg + newPinj - oldPinj
        refgen = zeros(size(ref));
        for k = 1:length(ref)
            temp = find(gbus == ref(k));
            refgen(k) = on(temp(1));
        end
        gen(refgen, PG) = gen(refgen, PG) + (B(ref, :) * Va - Pbus(ref)) * baseMVA;
    else                                %% AC formulation
        alg = upper(mpopt.pf.alg);
        if mpopt.verbose > 0
            switch alg
                case 'NR'
                    solver = 'Newton';
                case 'FDXB'
                    solver = 'fast-decoupled, XB';
                case 'FDBX'
                    solver = 'fast-decoupled, BX';
                case 'GS'
                    solver = 'Gauss-Seidel';
                case 'PQSUM'
                    solver = 'Power Summation';
                case 'ISUM'
                    solver = 'Current Summation';
                case 'YSUM'
                    solver = 'Admittance Summation';
                otherwise
                    solver = 'unknown';
            end
            fprintf(' -- AC Power Flow (%s)\n', solver);
        end
        %% initial state
        % V0    = ones(size(bus, 1), 1);            %% flat start
        V0  = bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));
        vcb = ones(size(V0));           %% create mask of voltage-controlled buses
        vcb(pq) = 0;                    %% exclude PQ buses
        k = find(vcb(gbus));            %% in-service gens at v-c buses
        V0(gbus(k)) = gen(on(k), VG) ./ abs(V0(gbus(k))).* V0(gbus(k));
    
        if qlim
            ref0 = ref;                         %% save index and angle of
            Varef0 = bus(ref0, VA);             %%   original reference bus(es)
            limited = [];                       %% list of indices of gens @ Q lims
            fixedQg = zeros(size(gen, 1), 1);   %% Qg of gens at Q limits
        end

        %% build admittance matrices
        [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
    
        repeat = 1;
        while (repeat)
            %% function for computing V dependent complex bus power injections
            %% (generation - load)
            Sbus = @(Vm)makeSbus(baseMVA, bus, gen, mpopt, Vm);
        
            %% run the power flow
            switch alg
                case 'NR'
                    [V, success, iterations] = newtonpf(Ybus, Sbus, V0, ref, pv, pq, mpopt);
                case {'FDXB', 'FDBX'}
                    [Bp, Bpp] = makeB(baseMVA, bus, branch, alg);
                    [V, success, iterations] = fdpf(Ybus, Sbus, V0, Bp, Bpp, ref, pv, pq, mpopt);
                case 'GS'
                    if (~isempty(mpopt.exp.sys_wide_zip_loads.pw) && ...
                            any(mpopt.exp.sys_wide_zip_loads.pw(2:3))) || ...
                            (~isempty(mpopt.exp.sys_wide_zip_loads.qw) && ...
                            any(mpopt.exp.sys_wide_zip_loads.qw(2:3)))
                        warning('runpf: Gauss-Seidel algorithm does not support ZIP load model. Converting to constant power loads.')
                        mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads', ...
                                        struct('pw', [], 'qw', []));
                    end
                    [V, success, iterations] = gausspf(Ybus, Sbus([]), V0, ref, pv, pq, mpopt);
                case {'PQSUM', 'ISUM', 'YSUM'}
                    [mpc, success, iterations] = radial_pf(mpc,mpopt);
                otherwise
                    error('runpf: ''%s'' is not a valid power flow algorithm. See ''pf.alg'' details in MPOPTION help.', alg);
            end
            its = its + iterations;
        
            %% update data matrices with solution
            switch alg
                case {'NR', 'FDXB', 'FDBX', 'GS'}
                    [bus, gen, branch] = pfsoln(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V, ref, pv, pq, mpopt);
                case {'PQSUM', 'ISUM', 'YSUM'}
                    bus = mpc.bus;
                    gen = mpc.gen;
                    branch = mpc.branch;
            end
        
            if success && qlim      %% enforce generator Q limits
                %% find gens with violated Q constraints
                mx = find( gen(:, GEN_STATUS) > 0 ...
                        & gen(:, QG) > gen(:, QMAX) + mpopt.opf.violation );
                mn = find( gen(:, GEN_STATUS) > 0 ...
                        & gen(:, QG) < gen(:, QMIN) - mpopt.opf.violation );
            
                if ~isempty(mx) || ~isempty(mn)  %% we have some Q limit violations
                    %% first check for INFEASIBILITY
                    infeas = union(mx', mn')';  %% transposes handle fact that
                        %% union of scalars is a row vector
                    remaining = find( gen(:, GEN_STATUS) > 0 & ...
                                    ( bus(gen(:, GEN_BUS), BUS_TYPE) == PV | ...
                                      bus(gen(:, GEN_BUS), BUS_TYPE) == REF ));
                    if length(infeas) == length(remaining) && all(infeas == remaining) && ...
                            (isempty(mx) || isempty(mn))
                        %% all remaining PV/REF gens are violating AND all are
                        %% violating same limit (all violating Qmin or all Qmax)
                        if mpopt.verbose
                            fprintf('All %d remaining gens exceed their Q limits : INFEASIBLE PROBLEM\n', length(infeas));
                        end
                        success = 0;
                        break;
                    end

                    %% one at a time?
                    if qlim == 2    %% fix largest violation, ignore the rest
                        [junk, k] = max([gen(mx, QG) - gen(mx, QMAX);
                                         gen(mn, QMIN) - gen(mn, QG)]);
                        if k > length(mx)
                            mn = mn(k-length(mx));
                            mx = [];
                        else
                            mx = mx(k);
                            mn = [];
                        end
                    end

                    if mpopt.verbose && ~isempty(mx)
                        fprintf('Gen %d at upper Q limit, converting to PQ bus\n', mx);
                    end
                    if mpopt.verbose && ~isempty(mn)
                        fprintf('Gen %d at lower Q limit, converting to PQ bus\n', mn);
                    end
                
                    %% save corresponding limit values
                    fixedQg(mx) = gen(mx, QMAX);
                    fixedQg(mn) = gen(mn, QMIN);
                    mx = [mx;mn];
                
                    %% convert to PQ bus
                    gen(mx, QG) = fixedQg(mx);      %% set Qg to binding limit
                    gen(mx, GEN_STATUS) = 0;        %% temporarily turn off gen,
                    for i = 1:length(mx)            %% (one at a time, since
                        bi = gen(mx(i), GEN_BUS);   %%  they may be at same bus)
                        bus(bi, [PD,QD]) = ...      %% adjust load accordingly,
                            bus(bi, [PD,QD]) - gen(mx(i), [PG,QG]);
                    end
                    if length(ref) > 1 && any(bus(gen(mx, GEN_BUS), BUS_TYPE) == REF)
                        error('runpf: Sorry, MATPOWER cannot enforce Q limits for slack buses in systems with multiple slacks.');
                    end
                    bus(gen(mx, GEN_BUS), BUS_TYPE) = PQ;   %% & set bus type to PQ
                
                    %% update bus index lists of each type of bus
                    ref_temp = ref;
                    [ref, pv, pq] = bustypes(bus, gen);
                    %% previous line can modify lists to select new REF bus
                    %% if there was none, so we should update bus with these
                    %% just to keep them consistent
                    if ref ~= ref_temp
                        bus(ref, BUS_TYPE) = REF;
                        bus( pv, BUS_TYPE) = PV;
                        if mpopt.verbose
                            fprintf('Bus %d is new slack bus\n', ref);
                        end
                    end
                    limited = [limited; mx];
                else
                    repeat = 0; %% no more generator Q limits violated
                end
            else
                repeat = 0;     %% don't enforce generator Q limits, once is enough
            end
        end
        if qlim && ~isempty(limited)
            %% restore injections from limited gens (those at Q limits)
            gen(limited, QG) = fixedQg(limited);    %% restore Qg value,
            for i = 1:length(limited)               %% (one at a time, since
                bi = gen(limited(i), GEN_BUS);      %%  they may be at same bus)
                bus(bi, [PD,QD]) = ...              %% re-adjust load,
                    bus(bi, [PD,QD]) + gen(limited(i), [PG,QG]);
            end
            gen(limited, GEN_STATUS) = 1;               %% and turn gen back on
            if ref ~= ref0
                %% adjust voltage angles to make original ref bus correct
                bus(:, VA) = bus(:, VA) - bus(ref0, VA) + Varef0;
            end
        end
    end
else
    t0 = tic;
    success = 0;
    its = 0;
    if mpopt.verbose
        fprintf('Power flow not valid : MATPOWER case contains no connected buses');
    end
end
mpc.et = toc(t0);
mpc.success = success;
mpc.iterations = its;

%%-----  output results  -----
%% convert back to original bus numbering & print results
[mpc.bus, mpc.gen, mpc.branch] = deal(bus, gen, branch);
results = int2ext(mpc);

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

if nargout == 1 || nargout == 2
    MVAbase = results;
    bus = success;
elseif nargout > 2
    [MVAbase, bus, gen, branch, et] = ...
        deal(results.baseMVA, results.bus, results.gen, results.branch, results.et);
% else  %% don't define MVAbase, so it doesn't print anything
end
