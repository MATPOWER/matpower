function [best_bus, best_gen, best_branch, best_f, best_success, et] = ...
		uopf(baseMVA, bus, gen, gencost, branch, area, mpopt)
%UOPF  Solves combined unit decommitment / optimal power flow.
%   [bus, gen, branch, f, success, et] = uopf(baseMVA, bus, gen, gencost, ...
%                   branch, area, mpopt)
%   Solves a combined unit decommitment and optimal power flow for a single
%   time period. Uses a heuristic that shuts off one generator at a time until
%   all generators dispatched at Pmin have Pmin * lambda greater than or equal
%   to the generator's cost of operating at Pmin (i.e. lambda >= avg cost).
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

%% do further decommitment as necessary
while 1
	count = count + 1;
	
	%% reassign bus types
	[ref, pv, pq] = bustypes(bus, gen);

	%% run opf
	%% turn down verbosity one level for call to opf
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
		
	if count == 1
		%% set "best" variables to new values
		best_bus		= bus;
		best_gen		= gen;
		best_branch		= branch;
		best_f			= f;
		best_success	= success;
		ignore			= find(gen(:, GEN_STATUS) <= 0);
		if ~success
			if verbose
				fprintf('UOPF: non-convergent OPF.\n');
			end
			break;
		end
	else
		if success
			if f < best_f
				%% made some progress, record solution
				best_bus		= bus;
				best_gen		= gen;
				best_branch		= branch;
				best_f			= f;
				best_success	= success;
				ignore			= find(gen(:, GEN_STATUS) <= 0);
			else
				%% return previous OPF solution
				if verbose
					fprintf('Cost did not decrease, restarting generator %d.\n', i);
				end

				%% restart generator i
				bus		= best_bus;
				gen		= best_gen;
				branch	= best_branch;
				f		= best_f;
				success	= best_success;
				ignore	= [ignore; i]; 	%% generator i no longer a candidate to be shut down
				di(ignore) = zeros(size(ignore));
	
				%% any other candidates to be shut down?
				if ~any(di > 1e-5)
					break;
				end
			end
		else
			%% OPF didn't converge, return previous commitment schedule
			if verbose
				fprintf('UOPF: non-convergent OPF, restarting generator %d.\n', i);
			end
			
			%% restart generator i
			bus		= best_bus;
			gen		= best_gen;
			branch	= best_branch;
			f		= best_f;
			success	= best_success;
			ignore	= [ignore; i]; 	%% generator i no longer a candidate to be shut down
			di(ignore) = zeros(size(ignore));

			%% any other candidates to be shut down?
			if ~any(di > 1e-5)
				break;
			end
		end
	end

	%% compute relative cost savings for shutting down each generator
	%% call it the "decommitment index"
	di = totcost(gencost, gen(:, PG)) - bus(gen(:, GEN_BUS), LAM_P) .* gen(:, PG);
	if sum(di > 0) > 4		%% compare prices alone, rather than price * qty
		ng = size(gen, 1);
		di = zeros(ng, 1);
		nzg = find(gen(:, PG));
		di(nzg) = totcost(gencost(nzg, :), gen(nzg, PG)) ./ gen(nzg, PG) - bus(gen(nzg, GEN_BUS), LAM_P);
	end
	di(ignore) = zeros(size(ignore));	%% generators which are no longer candidates
	if verbose > 1
		fprintf('\nDecommitment indices:\n');
		fprintf('\t%.4g\n', di);
	end
	
	%% check for units to decommit
	if any(di > 1e-5)
		%% shut down generator with largest decommitment index
		[junk, i] = fairmax(di);
		if verbose
			fprintf('Shutting down generator %d.\n', i);
		end

		%% set generation to zero
		gen(i, PG)			= 0;
		gen(i, QG)			= 0;
		gen(i, GEN_STATUS)	= 0;
	else
		break;
	end
end

%% compute elapsed time
et = etime(clock, t0);

return;
