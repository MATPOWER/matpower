function [MVAbase, bus, gen, gencost, branch, f, success, et] = ...
				runopf(casename, mpopt, fname)
%RUNOPF  Runs an optimal power flow.
%
%   [baseMVA, bus, gen, gencost, branch, f, success, et] = ...
%           runopf(casename, mpopt, fname)
%   
%   Runs an optimal power flow where casename is the name of the m-file
%   (without the .m extension) containing the opf data, and mpopt is a
%   MATPOWER options vector (see 'help mpoption' for details). Uses default
%   options if 2nd parameter is not given, and 'case' if 1st parameter
%   is not given. The results may optionally be printed to a file (appended
%   if the file exists) whose name is given in fname (in addition to
%   printing to STDOUT). Optionally returns the final values of baseMVA,
%   bus, gen, gencost, branch, f, success, and et.

%   MATPOWER Version 2.0
%   by Ray Zimmerman, PSERC Cornell    12/16/97
%   Copyright (c) 1996, 1997 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%% define named indices into data matrices
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

%% read data & convert to internal bus numbering
[baseMVA, bus, gen, branch, area, gencost] = feval(casename);
[i2e, bus, gen, branch, area] = ext2int(bus, gen, branch, area);

%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(bus, gen);

%% build admittance matrices
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% run the optimal power flow
[bus, gen, branch, f, success, et] = opf(baseMVA, bus, gen, gencost, branch, ...
					Ybus, Yf, Yt, ref, pv, pq, mpopt);

%% convert back to original bus numbering & print results
[bus, gen, branch, area] = int2ext(i2e, bus, gen, branch, area);
if fname
	[fd, msg] = fopen(fname, 'at');
	if fd == -1
		error(msg);
	else
		printpf(baseMVA, bus, gen, branch, f, success, et, fd, mpopt);
		fclose(fd);
	end
end
printpf(baseMVA, bus, gen, branch, f, success, et, 1, mpopt);

%% this is just to prevent it from printing baseMVA
%% when called with no output arguments
if nargout, MVAbase = baseMVA; end

return;
