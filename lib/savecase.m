function savecase(fname, p1, p2, p3, p4, p5, p6, p7)
%SAVECASE  Saves a MATPOWER case file, given a filename and the data matrices.
%
%   savecase(fname, baseMVA, bus, gen, branch, areas, gencost)
%       or
%   savecase(fname, comment, baseMVA, bus, gen, branch, areas, gencost)
%
%   Writes a MATPOWER case file, given a filename and the data matrices.
%   The fname parameter is the name of the file to be created or
%   overwritten. If fname ends with '.mat' it saves the case as a MAT-file
%   otherwise it saves it as an M-file.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman & Carlos Murillo, PSERC Cornell
%   Copyright (c) 1996-2003 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
	GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
	RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch;
[AREA_I, PRICE_REF_BUS] = idx_area;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, N, COST] = idx_cost;

%% default arguments
if isstr(p1)
  comment = p1;
  baseMVA = p2;
  bus     = p3;
  gen     = p4;
  branch  = p5;
  areas   = p6;
  gencost = p7;
else
  comment = '';
  baseMVA = p1;
  bus     = p2;
  gen     = p3;
  branch  = p4;
  areas   = p5;
  gencost = p6;
end

%% verify valid filename
l = length(fname);
rootname = [];
if l > 2
	if strcmp(fname(l-1:l), '.m')
		rootname = fname(1:l-2);
		extension = '.m';
	elseif l > 4
		if strcmp(fname(l-3:l), '.mat')
			rootname = fname(1:l-4);
			extension = '.mat';
		end
	end
end
if isempty(rootname)
	rootname = fname;
	extension = '.m';
	fname = [rootname, extension];
end

%% open and write the file
if strcmp(extension, '.mat')		%% MAT-file
	eval(['save ', rootname, ' baseMVA bus gen branch areas gencost;']);
else								%% M-file
	%% open file
	[fd, msg] = fopen(fname, 'wt');		%% print it to an m-file
	if fd == -1
		error(['savecase: ', msg]);
	end
	
	%% function header, etc.
	if isempty(areas) | isempty(gencost)
		fprintf(fd, 'function [baseMVA, bus, gen, branch] = %s\n\n', rootname);
	else
		fprintf(fd, 'function [baseMVA, bus, gen, branch, areas, gencost] = %s\n', rootname);
	end
	if length(comment) ~= 0
		fprintf(fd, '%% %s\n', comment);
	end
	fprintf(fd, '\n%%%%-----  Power Flow Data  -----%%%%\n');
	fprintf(fd, '%%%% system MVA base\n');
	fprintf(fd, 'baseMVA = %.4f;\n\n', baseMVA);
	
	%% bus data
	ncols = size(bus, 2);
	fprintf(fd, '%%%% bus data\n');
	fprintf(fd, '%%\ttype\tPd\tQd\tGs\tBs\tarea\tVm\tVa\tbaseKV\tzone\tVmax\tVmin\n');
	fprintf(fd, 'bus = [\n');
	if ncols < MU_VMIN				%% opf NOT SOLVED, save without lambda's & mu's
		fprintf(fd, '\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%.8f\t%.8f\t%.4f\t%d\t%.4f\t%.4f;\n', bus(:, 1:VMIN).');
	else							%% opf SOLVED, save with lambda's & mu's
		fprintf(fd, '\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%.8f\t%.8f\t%.4f\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f;\n', bus(:, 1:MU_VMIN).');
	end
	fprintf(fd, '];\n\n');
	
	%% generator data
	ncols = size(gen, 2);
	fprintf(fd, '%%%% generator data\n');
	fprintf(fd, '%%\tbus\tPg\tQg\tQmax\tQmin\tVg\tmBase\tstatus\tPmax\tPmin\n');
	fprintf(fd, 'gen = [\n');
	if ncols < MU_QMIN				%% opf NOT SOLVED, save without mu's
		fprintf(fd, '\t%d\t%.8f\t%.8f\t%.4f\t%.4f\t%.8f\t%.4f\t%d\t%.4f\t%.4f;\n', gen(:, 1:PMIN).');
	else							%% opf SOLVED, save with mu's
		fprintf(fd, '\t%d\t%.8f\t%.8f\t%.4f\t%.4f\t%.8f\t%.4f\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f;\n', gen(:, 1:MU_QMIN).');
	end
	fprintf(fd, '];\n\n');
	
	%% branch data
	ncols = size(branch, 2);
	fprintf(fd, '%%%% branch data\n');
	fprintf(fd, '%%\tfbus\ttbus\tr\tx\tb\trateA\trateB\trateC\tratio\tangle\tstatus\n');
	fprintf(fd, 'branch = [\n');
	if ncols < QT					%% power flow NOT SOLVED, save without line flows or mu's
		fprintf(fd, '\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%d;\n', branch(:, 1:BR_STATUS).');
	elseif ncols < MU_ST			%% power flow SOLVED, save with line flows but without mu's
		fprintf(fd, '\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%.4f\t%.4f\t%.4f\t%.4f;\n', branch(:, 1:QT).');
	else							%% opf SOLVED, save with lineflows & mu's
		fprintf(fd, '\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f;\n', branch(:, 1:MU_ST).');
	end
	fprintf(fd, '];\n\n');
	
	%% OPF data
	%% area data
	fprintf(fd, '%%%%-----  OPF Data  -----%%%%\n');
	fprintf(fd, '%%%% area data\n');
	fprintf(fd, 'areas = [\n');
	if ~isempty(areas)
		fprintf(fd, '\t%d\t%d;\n', areas(:, 1:PRICE_REF_BUS).');
	end
	fprintf(fd, '];\n\n');
		
	%% generator cost data
	fprintf(fd, '%%%% generator cost data\n');
	fprintf(fd, '%%\t1\tstartup\tshutdown\tn\tx0\ty0\t...\txn\tyn\n');
	fprintf(fd, '%%\t2\tstartup\tshutdown\tn\tc(n-1)\t...\tc0\n');
	fprintf(fd, 'gencost = [\n');
	if ~isempty(gencost)
		n = gencost(1, N);
		if gencost(1, MODEL) == PW_LINEAR
			n = 2 * n;
		end
		template = '\t%d\t%.2f\t%.2f\t%d';
		for i = 1:n
			template = [template, '\t%.6f'];
		end
		template = [template, ';\n'];
		fprintf(fd, template, gencost.');
	end
	fprintf(fd, '];\n\n');
	
	%% end
	fprintf(fd, 'return\n');
	
	%% close file
	if fd ~= 1
		fclose(fd);
	end
end

return;
