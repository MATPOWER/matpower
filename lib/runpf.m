function [MVAbase, bus, gen, branch, success, et] = runpf(casename, mpopt, fname, solvedcase)
%RUNPF  Runs a power flow.
%
%   [baseMVA, bus, gen, branch, success, et] = ...
%           runpf(casename, mpopt, fname, solvedcase)
%
%   Runs a power flow (full AC Newton's method by default) and optionally
%   returns the solved values in the data matrices, the objective function
%   value, a flag which is true if the algorithm was successful in finding a
%   solution, and the elapsed time in seconds. All input arguments are
%   optional. If casename is provided it specifies the name of the input
%   data file or struct (see also 'help caseformat' and 'help loadcase')
%   containing the power flow data. The default value is 'case9'. If the
%   mpopt is provided it overrides the default MATPOWER options vector and
%   can be used to specify the solution algorithm and output options among
%   other things (see 'help mpoption' for details). If the 3rd argument is
%   given the pretty printed output will be appended to the file whose name
%   is given in fname. If solvedcase is specified the solved case will be
%   written to a case file in MATPOWER format with the specified name. If
%   solvedcase ends with '.mat' it saves the case as a MAT-file otherwise it
%   saves it as an M-file.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2003 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%%-----  initialize  -----
%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
	RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
	GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;

%% default arguments
if nargin < 4
	solvedcase = '';				%% don't save solved case
	if nargin < 3
		fname = '';					%% don't print results to a file
		if nargin < 2
			mpopt = mpoption;		%% use default options
			if nargin < 1
				casename = 'case9';	%% default data file is 'case9.m'
			end
		end
	end
end

%% options
dc = mpopt(10);						%% use DC formulation?

%% read data & convert to internal bus numbering
[baseMVA, bus, gen, branch] = loadcase(casename);
[i2e, bus, gen, branch] = ext2int(bus, gen, branch);

%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(bus, gen);

%% generator info
on = find(gen(:, GEN_STATUS) > 0);		%% which generators are on?
gbus = gen(on, GEN_BUS);				%% what buses are they at?

%%-----  run the power flow  -----
t0 = clock;
if dc								%% DC formulation
	%% initial state
	Va0	= bus(:, VA) * (pi/180);
	
	%% build B matrices and phase shift injections
	[B, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);
	
	%% compute complex bus power injections (generation - load)
	%% adjusted for phase shifters and real shunts
	Pbus = real(makeSbus(baseMVA, bus, gen)) - Pbusinj - bus(:, GS) / baseMVA;
	
	%% "run" the power flow
	Va = dcpf(B, Pbus, Va0, ref, pv, pq);
	
	%% update data matrices with solution
	branch(:, [QF, QT]) = zeros(size(branch, 1), 2);
	branch(:, PF) = (Bf * Va + Pfinj) * baseMVA;
	branch(:, PT) = -branch(:, PF);
	bus(:, VM) = ones(size(bus, 1), 1);
	bus(:, VA) = Va * (180/pi);
	%% update Pg for swing generator (note: other gens at ref bus are accounted for in Pbus)
	%%		Pg = Pinj + Pload + Gs
	%%		newPg = oldPg + newPinj - oldPinj
	refgen = find(gbus == ref);				%% which is(are) the reference gen(s)?
	gen(on(refgen(1)), PG) = gen(on(refgen(1)), PG) + (B(ref, :) * Va - Pbus(ref)) * baseMVA;
	
	success = 1;
else								%% AC formulation
	%% initial state
	% V0	= ones(size(bus, 1), 1);			%% flat start
	V0	= bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));
	V0(gbus) = gen(on, VG) ./ abs(V0(gbus)).* V0(gbus);
	
	%% build admittance matrices
	[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
	
	%% compute complex bus power injections (generation - load)
	Sbus = makeSbus(baseMVA, bus, gen);
	
	%% run the power flow
	alg = mpopt(1);
	if alg == 1
		[V, success, iterations] = newtonpf(Ybus, Sbus, V0, ref, pv, pq, mpopt);
	elseif alg == 2 | alg == 3
		[Bp, Bpp] = makeB(baseMVA, bus, branch, alg);
		[V, success, iterations] = fdpf(Ybus, Sbus, V0, Bp, Bpp, ref, pv, pq, mpopt);
	elseif alg == 4
		[V, success, iterations] = gausspf(Ybus, Sbus, V0, ref, pv, pq, mpopt);
	else
		error('Only Newton''s method, fast-decoupled, and Gauss-Seidel power flow algorithms currently implemented.');
	end
	
	%% update data matrices with solution
	[bus, gen, branch] = pfsoln(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V, ref, pv, pq);
end
et = etime(clock, t0);

%%-----  output results  -----
%% convert back to original bus numbering & print results
[bus, gen, branch] = int2ext(i2e, bus, gen, branch);
if fname
	[fd, msg] = fopen(fname, 'at');
	if fd == -1
		error(msg);
	else
		printpf(baseMVA, bus, gen, branch, [], success, et, fd, mpopt);
		fclose(fd);
	end
end
printpf(baseMVA, bus, gen, branch, [], success, et, 1, mpopt);

%% save solved case
if solvedcase
	savecase(solvedcase, baseMVA, bus, gen, branch);
end

%% this is just to prevent it from printing baseMVA
%% when called with no output arguments
if nargout, MVAbase = baseMVA; end

return;
