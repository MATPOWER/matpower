function [best_bus, best_gen, best_branch, best_f, best_success, et] = ...
		uopf(baseMVA, bus, gen, gencost, branch, Ybus, Yf, Yt, mpopt)
%UOPF  Solves combined unit commitment / optimal power flow.
%   [bus, gen, branch, f, success, et] = uopf(baseMVA, bus, gen, gencost, ...
%                   branch, Ybus, Yf, Yt, mpopt)
%   Assumes that Ybus, Yf, Yt are consistent with (or override) data in bus,
%   gen, branch. If VERBOSE in mpopt (see 'help mpoption') is true, it prints
%   progress info, if it's > 1 it prints the output of each individual opf.

%   MATPOWER Version 2.0
%   by Ray Zimmerman, PSERC Cornell    12/11/97
%   Copyright (c) 1996, 1997 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%%----- initialization -----
count		= 0;
i			= 0;		%% this is to work around a bug in Matlab (4 and 5)

%% default arguments
if nargin < 9
	mpopt = mpoption;		%% use default options
end

%% options
verbose	= mpopt(31);

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
	GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, N, COST] = idx_cost;

%%-----  do combined unit commitment/optimal power flow  -----
t0 = clock;									%% start timer

%% do decommitment as necessary
count = 0;
while 1
	count = count + 1;
	
	%% reassign bus types
	[ref, pv, pq] = bustypes(bus, gen);

	%% run opf
	%% turn down verbosity one level for call to opf
	if verbose
		mpopt = mpoption(mpopt, 'VERBOSE', verbose-1);
	end
	[bus, gen, branch, f, success, et] =  opf(baseMVA, bus, gen, gencost, branch, ...
					Ybus, Yf, Yt, ref, pv, pq, mpopt);
		
	if count == 1
		%% set "best" variables to new values
		best_bus		= bus;
		best_gen		= gen;
		best_branch		= branch;
		best_f			= f;
		best_success	= success;
		ignore			= find(gen(:, GEN_STATUS) == 0);
		if ~success
			if verbose
				fprintf('\nUOFP: non-convergent OPF.\n');
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
				ignore			= find(gen(:, GEN_STATUS) == 0);
			else
				%% return previous OPF solution
				if verbose
					fprintf('\nCost did not decrease, restarting generator %d.\n', i);
				end

				%% restart generator i
				bus		= best_bus;
				gen		= best_gen;
				branch	= best_branch;
				f		= best_f;
				success	= best_success;
				ignore	= find(gen(:, GEN_STATUS) == 0);
				ignore	= [ignore, i]; 	%% generator i no longer a candidate to be shut down
				di(ignore) = zeros(size(ignore));
	
				%% any other candidates to be shut down?
				if ~any(di > 1e-5)
					break;
				end
			end
		else
			%% OPF didn't converge, return previous commitment schedule
			if verbose
				fprintf('\nUOFP: non-convergent OPF, restarting generator %d.\n', i);
			end
			
			%% restart generator i
			bus		= best_bus;
			gen		= best_gen;
			branch	= best_branch;
			f		= best_f;
			success	= best_success;
			ignore	= find(gen(:, GEN_STATUS) == 0);
			ignore	= [ignore, i]; 	%% generator i no longer a candidate to be shut down
			di(ignore) = zeros(size(ignore));

			%% any other candidates to be shut down?
			if ~any(di > 1e-5)
				break;
			end
		end
	end

	%% compute relative cost savings for shutting down each generator
	%% call it the "decommitment index"
	di = totcost(gencost, gen(:, PG)) -  bus(gen(:, GEN_BUS), LAM_P) .* gen(:, PG);
	di(ignore) = zeros(size(ignore));	%% generators which are no longer candidates
	if verbose > 1
		fprintf('\nDecommitment indices:\n');
		fprintf('\t%.4g\n', di);
	end
	
	%% check for units to decommit
	if any(di > 1e-5)
		%% shut down generator with largest decommitment index
		[junk, i] = max(di);
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
