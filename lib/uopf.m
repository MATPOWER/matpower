function [bus, gen, branch, f, success, info, et, g, jac, xr, pimul] = ...
    uopf(varargin)
%UOPF  Solves combined unit decommitment / optimal power flow.
%
%   Returns either a results struct and optinally a success flag, or individual
%   data matrices, the objective function value and a success flag. In the
%   latter case, there are additional optional return values.
%
%   results = uopf(...)
%   [results, success] = uopf(...)
%   [bus, gen, branch, f, success] = uopf(...)
%   [bus, gen, branch, f, success, info, et, g, jac, xr, pimul] = uopf(...)
%
%   Input arguments options are as follows:
%
%   uopf(mpc)
%   uopf(mpc, mpopt)
%   uopf(mpc, userfcn, mpopt)
%   uopf(mpc, A, l, u)
%   uopf(mpc, A, l, u, mpopt)
%   uopf(mpc, A, l, u, mpopt, N, fparm, H, Cw)
%   uopf(mpc, A, l, u, mpopt, N, fparm, H, Cw, z0, zl, zu)
%
%   uopf(baseMVA, bus, gen, branch, areas, gencost)
%   uopf(baseMVA, bus, gen, branch, areas, gencost, mpopt)
%   uopf(baseMVA, bus, gen, branch, areas, gencost, userfcn, mpopt)
%   uopf(baseMVA, bus, gen, branch, areas, gencost, A, l, u)
%   uopf(baseMVA, bus, gen, branch, areas, gencost, A, l, u, mpopt)
%   uopf(baseMVA, bus, gen, branch, areas, gencost, A, l, u, ...
%                               mpopt, N, fparm, H, Cw)
%   uopf(baseMVA, bus, gen, branch, areas, gencost, A, l, u, ...
%                               mpopt, N, fparm, H, Cw, z0, zl, zu)
%
%   See 'help opf' for more information on input and output arguments.
%
%   Solves a combined unit decommitment and optimal power flow for a single
%   time period. Uses an algorithm similar to dynamic programming. It proceeds
%   through a sequence of stages, where stage N has N generators shut down,
%   starting with N=0. In each stage, it forms a list of candidates (gens at
%   their Pmin limits) and computes the cost with each one of them shut down.
%   It selects the least cost case as the starting point for the next stage,
%   continuing until there are no more candidates to be shut down or no
%   more improvement can be gained by shutting something down.
%   If VERBOSE in mpopt (see 'help mpoption') is true, it prints progress
%   info, if it is > 1 it prints the output of each individual opf.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialization -----
t0 = clock;                                 %% start timer

%% process input arguments
[mpc, mpopt] = opf_args(varargin{:});

%% options
verbose = mpopt(31);    %% VERBOSE
if verbose      %% turn down verbosity one level for calls to opf
    mpopt = mpoption(mpopt, 'VERBOSE', verbose-1);
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

    if verbose
        fprintf('Shutting down generator %d so all Pmin limits can be satisfied.\n', i);
    end

    %% set generation to zero
    mpc.gen(i, [ PG QG GEN_STATUS ]) = 0;
    
    %% update minimum gen capacity
    on  = find( mpc.gen(:, GEN_STATUS) > 0 & ~isload(mpc.gen) );   %% gens in service
    Pmin = mpc.gen(on, PMIN);
end

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
        [results, success] = opf(mpc, mpopt);
        
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
        if verbose
            fprintf('Shutting down generator %d.\n', k1);
        end
        
        results0 = results1;
        mpc.bus = results0.bus;     %% use these V as starting point for OPF
    end 
end

%% compute elapsed time
et = etime(clock, t0);

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
  printpf(results0, et, 1, mpopt);
end
