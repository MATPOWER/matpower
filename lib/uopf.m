function [bus0, gen0, branch0, f0, success0, et] = ...
		uopf(baseMVA, bus, gen, gencost, branch, area, mpopt)
%UOPF  Solves combined unit decommitment / optimal power flow.
%   [bus, gen, branch, f, success, et] = uopf(baseMVA, bus, gen, gencost, ...
%                   branch, area, mpopt)
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
%   Copyright (c) 1996-2003 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%%----- initialization -----
count		= 0;
i			= 0;	%% this is to work around a bug in Matlab (4 and 5)

%% default arguments
if nargin < 6
	mpopt = mpoption;					%% use default options
end

%% options
verbose	= mpopt(31);
dc = mpopt(10);							%% use DC formulation?

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
	GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, N, COST] = idx_cost;

%%-----  do combined unit commitment/optimal power flow  -----
t0 = clock;									%% start timer

%% build network matrices
if dc								%% DC formulation
	%% build B matrices and phase shift injections
	[Bbus, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);
else								%% AC formulation
	%% build admittance matrices
	[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
end


%% check for sum(Pmin) > total load, decommit as necessary
load_capacity	= sum(bus(:, PD));		%% compute total load capacity
on = find(gen(:, GEN_STATUS) > 0);		%% which generators are on?
Pmin = gen(on, PMIN);
while sum(Pmin) > load_capacity
	%% shut down most expensive unit
	avgPmincost = totcost(gencost(on, :), Pmin) ./ Pmin;
	[junk, i] = fairmax(avgPmincost);	%% pick one with max avg cost at Pmin
	i = on(i);							%% convert to generator index

	if verbose
		fprintf('Shutting down generator %d so all Pmin limits can be satisfied.\n', i);
	end

	%% set generation to zero
	gen(i, PG)			= 0;
	gen(i, QG)			= 0;
	gen(i, GEN_STATUS)	= 0;
	
	%% update minimum gen capacity
	on = find(gen(:, GEN_STATUS) > 0);		%% which generators are on?
	Pmin = gen(on, PMIN);
end

%% run initial opf
%% turn down verbosity one level for call to opf
[ref, pv, pq] = bustypes(bus, gen);
if verbose
	mpopt = mpoption(mpopt, 'VERBOSE', verbose-1);
end
if dc								%% DC formulation
	[bus, gen, branch, f, success, et] = dcopf(baseMVA, bus, gen, gencost, branch, ...
						Bbus, Bf, Pbusinj, Pfinj, ref, pv, pq, mpopt);
else								%% AC formulation
	[bus, gen, branch, f, success, et] = opf(baseMVA, bus, gen, gencost, branch, ...
						area, Ybus, Yf, Yt, ref, pv, pq, mpopt);
end

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
	done = 1;	%% do not check for further decommitment unless we
				%%  see something better during this stage
	for i = 1:length(candidates)
		k = candidates(i);
		%% start with best for this stage
		gen = gen0;
		
		%% shut down gen k
		gen(k, PG)			= 0;
		gen(k, QG)			= 0;
		gen(k, GEN_STATUS)	= 0;
		
		%% run opf
		%% turn down verbosity one level for call to opf
		[ref, pv, pq] = bustypes(bus, gen);
		if verbose
			mpopt = mpoption(mpopt, 'VERBOSE', verbose-1);
		end
		if dc								%% DC formulation
			[bus, gen, branch, f, success, et] = dcopf(baseMVA, bus0, gen, gencost, branch0, ...
								Bbus, Bf, Pbusinj, Pfinj, ref, pv, pq, mpopt);
		else								%% AC formulation
			[bus, gen, branch, f, success, et] = opf(baseMVA, bus0, gen, gencost, branch0, ...
								area, Ybus, Yf, Yt, ref, pv, pq, mpopt);
		end
		
		%% something better?
		if success & f < f1
			bus1 = bus;
			gen1 = gen;
			branch1 = branch;
			success1 = success;
			f1 = f;
			k1 = k;
			done = 0;	%% make sure we check for further decommitment
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
