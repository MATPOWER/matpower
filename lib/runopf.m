function [MVAbase, bus, gen, gencost, branch, f, success, et] = ...
				runopf(casename, mpopt, fname, solvedcase)
%RUNOPF  Runs an optimal power flow.
%
%   [baseMVA, bus, gen, gencost, branch, f, success, et] = ...
%           runopf(casename, mpopt, fname, solvedcase)
%   
%   Runs an optimal power flow where casename is the name of the m-file
%   (without the .m extension) containing the opf data, and mpopt is a
%   MATPOWER options vector (see 'help mpoption' for details). Uses default
%   options if 2nd parameter is not given, and 'case9' if 1st parameter
%   is not given. The results may optionally be printed to a file (appended
%   if the file exists) whose name is given in fname (in addition to
%   printing to STDOUT). Optionally returns the final values of baseMVA,
%   bus, gen, gencost, branch, f, success, and et. If a name is given in
%   solvedcase, the solved case will be written to a case file in MATPOWER
%   format with the specified name with a '.m' extension added.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2003 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%%-----  initialize  -----
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
[baseMVA, bus, gen, branch, areas, gencost] = loadcase(casename);
[i2e, bus, gen, branch, areas] = ext2int(bus, gen, branch, areas);

%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(bus, gen);

%%-----  run the optimal power flow  -----
if dc								%% DC formulation
	%% build B matrices and phase shift injections
	[Bbus, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);
	
	%% run the optimal power flow
	[bus, gen, branch, f, success, et] = dcopf(baseMVA, bus, gen, gencost, branch, ...
						Bbus, Bf, Pbusinj, Pfinj, ref, pv, pq, mpopt);
else								%% AC formulation
	%% build admittance matrices
	[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
	
	%% run the optimal power flow
	[bus, gen, branch, f, success, et] = opf(baseMVA, bus, gen, gencost, branch, ...
						areas, Ybus, Yf, Yt, ref, pv, pq, mpopt);
end

%%-----  output results  -----
%% convert back to original bus numbering & print results
[bus, gen, branch, areas] = int2ext(i2e, bus, gen, branch, areas);
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

%% save solved case
if solvedcase
	savecase(solvedcase, baseMVA, bus, gen, branch, areas, gencost);
end

%% this is just to prevent it from printing baseMVA
%% when called with no output arguments
if nargout, MVAbase = baseMVA; end

return;
