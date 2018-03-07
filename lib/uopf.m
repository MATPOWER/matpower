function [bus, gen, branch, f, success, info, et, g, jac, xr, pimul] = ...
    uopf(varargin)
%UOPF  Solves combined unit decommitment / optimal power flow.
%   [RESULTS, SUCCESS] = UOPF(MPC, MPOPT)
%
%   Returns either a RESULTS struct and an optional SUCCESS flag, or individual
%   data matrices, the objective function value and a SUCCESS flag. In the
%   latter case, there are additional optional return values. See Examples
%   below for the possible calling syntax options.
%
%   Examples:
%       Output argument options:
%
%       results = uopf(...)
%       [results, success] = uopf(...)
%       [bus, gen, branch, f, success] = uopf(...)
%       [bus, gen, branch, f, success, info, et, g, jac, xr, pimul] = uopf(...)
%
%       Input arguments options:
%
%       uopf(mpc)
%       uopf(mpc, mpopt)
%       uopf(mpc, userfcn, mpopt)
%       uopf(mpc, A, l, u)
%       uopf(mpc, A, l, u, mpopt)
%       uopf(mpc, A, l, u, mpopt, N, fparm, H, Cw)
%       uopf(mpc, A, l, u, mpopt, N, fparm, H, Cw, z0, zl, zu)
%
%       uopf(baseMVA, bus, gen, branch, areas, gencost)
%       uopf(baseMVA, bus, gen, branch, areas, gencost, mpopt)
%       uopf(baseMVA, bus, gen, branch, areas, gencost, userfcn, mpopt)
%       uopf(baseMVA, bus, gen, branch, areas, gencost, A, l, u)
%       uopf(baseMVA, bus, gen, branch, areas, gencost, A, l, u, mpopt)
%       uopf(baseMVA, bus, gen, branch, areas, gencost, A, l, u, ...
%                                   mpopt, N, fparm, H, Cw)
%       uopf(baseMVA, bus, gen, branch, areas, gencost, A, l, u, ...
%                                   mpopt, N, fparm, H, Cw, z0, zl, zu)
%
%   See OPF for more information on input and output arguments.
%
%   Solves a combined unit decommitment and optimal power flow for a single
%   time period. Uses an algorithm similar to dynamic programming. It proceeds
%   through a sequence of stages, where stage N has N generators shut down,
%   starting with N=0. In each stage, it forms a list of candidates (gens at
%   their Pmin limits) and computes the cost with each one of them shut down.
%   It selects the least cost case as the starting point for the next stage,
%   continuing until there are no more candidates to be shut down or no
%   more improvement can be gained by shutting something down.
%   If MPOPT.verbose (see MPOPTION) is true, it prints progress
%   info, if it is > 1 it prints the output of each individual opf.
%
%   See also OPF, RUNUOPF.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialization -----
t0 = tic;       %% start timer

%% process input arguments
[mpc, mpopt] = opf_args(varargin{:});

%% options
if mpopt.verbose    %% turn down verbosity one level for calls to opf
    mpopt = mpoption(mpopt, 'verbose', mpopt.verbose-1);
end

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%%-----  do combined unit commitment/optimal power flow  -----

%% check for sum(Pmin) > total load, decommit as necessary
on   = find( mpc.gen(:, GEN_STATUS) > 0 & ~isload(mpc.gen) );   %% gens in service
onld = find( mpc.gen(:, GEN_STATUS) > 0 &  isload(mpc.gen) );   %% disp loads in serv
load_capacity = sum(mpc.bus(:, PD)) - sum(mpc.gen(onld, PMIN)); %% total load capacity
Pmin = mpc.gen(on, PMIN);
while sum(Pmin) > load_capacity
    %% shut down most expensive unit
    avgPmincost = totcost(mpc.gencost(on, :), Pmin) ./ Pmin;
    [junk, i] = fairmax(avgPmincost);   %% pick one with max avg cost at Pmin
    i = on(i);                          %% convert to generator index

    if mpopt.verbose
        fprintf('Shutting down generator %d so all Pmin limits can be satisfied.\n', i);
    end

    %% set generation to zero
    mpc.gen(i, [ PG QG GEN_STATUS ]) = 0;
    
    %% update minimum gen capacity
    on  = find( mpc.gen(:, GEN_STATUS) > 0 & ~isload(mpc.gen) );   %% gens in service
    Pmin = mpc.gen(on, PMIN);
end
if ~any(mpc.gen(:, GEN_STATUS) > 0)     %% don't bother to run anything if
    success = 0;                        %% everything has been shut down
    results0 = mpc;
    results0.success = success;
    results0.f = NaN;
    results0.et = 0;
    if mpopt.verbose
        fprintf('Infeasible problem, Pmin limits cannot be satisfied without shutting down all generators.\n');
    end
else
    %% run initial opf
    [results, success] = opf(mpc, mpopt);

    %% best case so far
    results1 = results;

    %% best case for this stage (ie. with n gens shut down, n=0,1,2 ...)
    results0 = results1;
    mpc.bus = results0.bus;     %% use these V as starting point for OPF

    while 1
        %% get candidates for shutdown
        candidates = find(results0.gen(:, MU_PMIN) > 0 & results0.gen(:, PMIN) > 0);
        if isempty(candidates)
            break;
        end
        done = 1;   %% do not check for further decommitment unless we
                    %%  see something better during this stage
        for i = 1:length(candidates)
            k = candidates(i);
            %% start with best for this stage
            mpc.gen = results0.gen;
        
            %% shut down gen k
            mpc.gen(k, [ PG QG GEN_STATUS ]) = 0;
        
            %% run opf
            if any(mpc.gen(:, GEN_STATUS) > 0)
                [results, success] = opf(mpc, mpopt);
            else
                success = 0;
            end
        
            %% something better?
            if success && results.f < results1.f
                results1 = results;
                k1 = k;
                done = 0;   %% make sure we check for further decommitment
            end
        end

        if done
            %% decommits at this stage did not help, so let's quit
            break;
        else
            %% shutting something else down helps, so let's keep going
            if mpopt.verbose
                fprintf('Shutting down generator %d.\n', k1);
            end
        
            results0 = results1;
            mpc.bus = results0.bus;     %% use these V as starting point for OPF
        end 
    end
end

%% compute elapsed time
et = toc(t0);

%% finish preparing output
if nargout > 0
  success = results0.success;
  if nargout <= 2
    results0.et = et;
    bus = results0;
    gen = success;
  else
    [bus, gen, branch, f, info, xr, pimul] = deal(results0.bus, results0.gen, ...
                    results0.branch, results0.f, results0.raw.info, ...
                    results0.raw.xr, results0.raw.pimul);
    if isfield(results0, 'g')
      g = results0.g;
    end
    if isfield(results0, 'dg')
      jac = results0.dg;
    end
  end
elseif results0.success
  results0.et = et;
  printpf(results0, 1, mpopt);
end


function [val, idx] = fairmax(x)
%FAIRMAX    Same as built-in MAX, except breaks ties randomly.
%   [VAL, IDX] = FAIRMAX(X) takes a vector as an argument and returns
%   the same output as the built-in function MAX with two output
%   parameters, except that where the maximum value occurs at more
%   than one position in the  vector, the index is chosen randomly
%   from these positions as opposed to just choosing the first occurance.
%
%   See also MAX.

val = max(x);               %% find max value
i   = find(x == val);       %% find all positions where this occurs
n   = length(i);            %% number of occurences
idx = i( fix(n*rand)+1 );   %% select index randomly among occurances
