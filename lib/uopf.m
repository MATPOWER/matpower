function [bus0, gen0, branch0, f0, success0, et] = ...
        uopf(baseMVA, bus, gen, gencost, branch, areas, mpopt)
%UOPF  Solves combined unit decommitment / optimal power flow.
%   [bus, gen, branch, f, success, et] = uopf(baseMVA, bus, gen, gencost, ...
%                   branch, areas, mpopt)
%   Solves a combined unit decommitment and optimal power flow for a single
%   time period. Uses an algorithm similar to dynamic programming. It proceeds
%   through a sequence of stages, where stage N has N generators shut down,
%   starting with N=0. In each stage, it forms a list of candidates (gens at
%   their Pmin limits) and computes the cost with each one of them shut down.
%   It selects the least cost case as the starting point for the next stage,
%   continuing until there are no more candidates to be shut down or no
%   more improvement can be gained by shutting something down.
%   If VERBOSE in mpopt (see 'help mpoption') is true, it prints progress
%   info, if it's > 1 it prints the output of each individual opf.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialization -----
count       = 0;
i           = 0;    %% this is to work around a bug in Matlab (4 and 5)

%% default arguments
if nargin < 6
    mpopt = mpoption;                   %% use default options
end

%% options
verbose = mpopt(31);
if verbose      %% turn down verbosity one level for calls to opf
    mpopt = mpoption(mpopt, 'VERBOSE', verbose-1);
end

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, ...
    PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, QMAX2, QMIN2, ...
    RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, N, COST] = idx_cost;

%%-----  do combined unit commitment/optimal power flow  -----
t0 = clock;                                 %% start timer

%% check for sum(Pmin) > total load, decommit as necessary
on   = find( gen(:, GEN_STATUS) > 0 & ~isload(gen) );   %% gens in service
onld = find( gen(:, GEN_STATUS) > 0 &  isload(gen) );   %% disp loads in serv
load_capacity = sum(bus(:, PD)) - sum(gen(onld, PMIN)); %% total load capacity
Pmin = gen(on, PMIN);
while sum(Pmin) > load_capacity
    %% shut down most expensive unit
    avgPmincost = totcost(gencost(on, :), Pmin) ./ Pmin;
    [junk, i] = fairmax(avgPmincost);   %% pick one with max avg cost at Pmin
    i = on(i);                          %% convert to generator index

    if verbose
        fprintf('Shutting down generator %d so all Pmin limits can be satisfied.\n', i);
    end

    %% set generation to zero
    gen(i, PG)          = 0;
    gen(i, QG)          = 0;
    gen(i, GEN_STATUS)  = 0;
    
    %% update minimum gen capacity
    on  = find( gen(:, GEN_STATUS) > 0 & ~isload(gen) );   %% gens in service
    Pmin = gen(on, PMIN);
end

%% run initial opf
[bus, gen, branch, f, success, info, et] = opf(baseMVA, bus, gen, branch, ...
                                   areas, gencost, mpopt);

%% best case so far
bus1 = bus;
gen1 = gen;
branch1 = branch;
success1 = success;
f1 = f;

%% best case for this stage (ie. with n gens shut down, n=0,1,2 ...)
bus0 = bus1;
gen0 = gen1;
branch0 = branch1;
success0 = success1;
f0 = f1;

while 1
    %% get candidates for shutdown
    candidates = find(gen0(:, MU_PMIN) > 0 & gen0(:, PMIN) > 0);
    if isempty(candidates)
        break;
    end
    done = 1;   %% do not check for further decommitment unless we
                %%  see something better during this stage
    for i = 1:length(candidates)
        k = candidates(i);
        %% start with best for this stage
        gen = gen0;
        
        %% shut down gen k
        gen(k, PG)          = 0;
        gen(k, QG)          = 0;
        gen(k, GEN_STATUS)  = 0;
        
        %% run opf
        [bus, gen, branch, f, success, info, et] = opf(baseMVA, bus0, gen, branch0, ...
                                           areas, gencost, mpopt);
        
        %% something better?
        if success & f < f1
            bus1 = bus;
            gen1 = gen;
            branch1 = branch;
            success1 = success;
            f1 = f;
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
        
        bus0 = bus1;
        gen0 = gen1;
        branch0 = branch1;
        success0 = success1;
        f0 = f1;
    end 
end

%% compute elapsed time
et = etime(clock, t0);

return;
