function [MVAbase, bus, gen, branch, success, et] = runpf(casename, mpopt, fname)
%RUNPF  Runs a power flow.
%
%   [baseMVA, bus, gen, branch, success, et] = runpf(casename, mpopt, fname)
%
%   Runs a full Newton's method power flow where casename is the name of
%   the m-file (without the .m extension) containing the power flow data,
%   and mpopt is a MATPOWER options vector (see 'help mpoption' for details).
%   Uses default options if 2nd parameter is not given, and 'case' if 1st
%   parameter is not given. The results may optionally be printed to a file
%   (appended if the file exists) whose name is given in fname (in addition
%   to printing to STDOUT). Optionally returns the final values of baseMVA,
%   bus, gen, branch, success, and et.

%   MATPOWER Version 2.0
%   by Ray Zimmerman, PSERC Cornell    12/24/97
%   Copyright (c) 1996, 1997 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.
tic;
%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
	GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;

%% default arguments
if nargin < 3
	fname = '';					%% don't print results to a file
	if nargin < 2
		mpopt = mpoption;		%% use default options
		if nargin < 1
			casename = 'case';	%% default data file is 'case.m'
		end
	end
end

%% options
alg = mpopt(1);

%% read data & convert to internal bus numbering
[baseMVA, bus, gen, branch] = feval(casename);
[i2e, bus, gen, branch] = ext2int(bus, gen, branch);

%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(bus, gen);

%% build admittance matrices
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% compute complex bus power injections (generation - load)
Sbus = makeSbus(baseMVA, bus, gen);

%% generator info
on = find(gen(:, GEN_STATUS));				%% which generators are on?

%% initialize V and Pg from data from case file
% V0	= ones(size(bus, 1), 1);		%% flat start
V0	= bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));
V0(gen(on, GEN_BUS)) = gen(on, VG) ./ abs(V0(gen(on, GEN_BUS))).* V0(gen(on, GEN_BUS));

%% run the power flow
t0 = clock;
if alg == 1
	[V, success, iterations] = newtonpf(Ybus, Sbus, V0, ref, pv, pq, mpopt);
elseif alg == 2 | alg == 3
	[Bp, Bpp] = makeB(baseMVA, bus, branch, alg);
	[V, success, iterations] = fdpf(Ybus, Sbus, V0, Bp, Bpp, ref, pv, pq, mpopt);
else
	error('Only Newton''s method and fast-decoupled power flow algorithms currently implemented.');
end

%% compute flows etc.
[bus, gen, branch] = pfsoln(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V, ref, pv, pq);
et = etime(clock, t0);

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

%% this is just to prevent it from printing baseMVA
%% when called with no output arguments
if nargout, MVAbase = baseMVA; end

return;
