function t_loadcase(quiet)
%T_LOADCASE  Test that loadcase() works with a struct as well as case file.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

if nargin < 1
	quiet = 0;
end

t_begin(53, quiet);

%% compare result of loading from m-file file to result of using data matrices
casefile = 't_case9_opf';
matfile  = 't_mat9_opf';
pfcasefile = 't_case9_pf';
pfmatfile  = 't_mat9_pf';
[baseMVA, bus, gen, branch, areas, gencost] = feval(casefile);

%% save as .mat file
eval(['save ' matfile ' baseMVA bus gen branch areas gencost']);
eval(['save ' pfmatfile ' baseMVA bus gen branch']);

t = 'loadcase(opf_M_file) without .m extension : ';
[baseMVA1, bus1, gen1, branch1, area1, gencost1] = loadcase(casefile);
t_is(baseMVA1,	baseMVA,	12, [t 'baseMVA']);
t_is(bus1,		bus,		12, [t 'bus']);
t_is(gen1,		gen,		12, [t 'gen']);
t_is(branch1,	branch,		12, [t 'branch']);
t_is(area1,		areas,		12, [t 'areas']);
t_is(gencost1,	gencost,	12, [t 'gencost']);

t = 'loadcase(opf_M_file) with .m extension : ';
[baseMVA1, bus1, gen1, branch1, area1, gencost1] = loadcase([casefile '.m']);
t_is(baseMVA1,	baseMVA,	12, [t 'baseMVA']);
t_is(bus1,		bus,		12, [t 'bus']);
t_is(gen1,		gen,		12, [t 'gen']);
t_is(branch1,	branch,		12, [t 'branch']);
t_is(area1,		areas,		12, [t 'areas']);
t_is(gencost1,	gencost,	12, [t 'gencost']);

t = 'load_case(opf_MAT_file) without .mat extension : ';
[baseMVA1, bus1, gen1, branch1, area1, gencost1] = loadcase(matfile);
t_is(baseMVA1,	baseMVA,	12, [t 'baseMVA']);
t_is(bus1,		bus,		12, [t 'bus']);
t_is(gen1,		gen,		12, [t 'gen']);
t_is(branch1,	branch,		12, [t 'branch']);
t_is(area1,		areas,		12, [t 'areas']);
t_is(gencost1,	gencost,	12, [t 'gencost']);

t = 'load_case(opf_MAT_file) with .mat extension : ';
[baseMVA1, bus1, gen1, branch1, area1, gencost1] = loadcase([matfile '.mat']);
t_is(baseMVA1,	baseMVA,	12, [t 'baseMVA']);
t_is(bus1,		bus,		12, [t 'bus']);
t_is(gen1,		gen,		12, [t 'gen']);
t_is(branch1,	branch,		12, [t 'branch']);
t_is(area1,		areas,		12, [t 'areas']);
t_is(gencost1,	gencost,	12, [t 'gencost']);

t = 'loadcase(pf_M_file) without .m extension : ';
[baseMVA1, bus1, gen1, branch1] = loadcase(pfcasefile);
t_is(baseMVA1,	baseMVA,	12, [t 'baseMVA']);
t_is(bus1,		bus,		12, [t 'bus']);
t_is(gen1,		gen,		12, [t 'gen']);
t_is(branch1,	branch,		12, [t 'branch']);

t = 'loadcase(pf_M_file) with .m extension : ';
[baseMVA1, bus1, gen1, branch1] = loadcase([pfcasefile '.m']);
t_is(baseMVA1,	baseMVA,	12, [t 'baseMVA']);
t_is(bus1,		bus,		12, [t 'bus']);
t_is(gen1,		gen,		12, [t 'gen']);
t_is(branch1,	branch,		12, [t 'branch']);

t = 'load_case(pf_MAT_file) without .mat extension : ';
[baseMVA1, bus1, gen1, branch1] = loadcase(pfmatfile);
t_is(baseMVA1,	baseMVA,	12, [t 'baseMVA']);
t_is(bus1,		bus,		12, [t 'bus']);
t_is(gen1,		gen,		12, [t 'gen']);
t_is(branch1,	branch,		12, [t 'branch']);

t = 'load_case(pf_MAT_file) with .mat extension : ';
[baseMVA1, bus1, gen1, branch1] = loadcase([pfmatfile '.mat']);
t_is(baseMVA1,	baseMVA,	12, [t 'baseMVA']);
t_is(bus1,		bus,		12, [t 'bus']);
t_is(gen1,		gen,		12, [t 'gen']);
t_is(branch1,	branch,		12, [t 'branch']);

delete([ matfile '.mat' ]);
delete([ pfmatfile '.mat' ]);

t = 'loadcase(my_struct) : ';
c.baseMVA	= baseMVA;
c.bus		= bus;
c.gen		= gen;
c.branch	= branch;
c.areas		= areas;
c.gencost	= gencost;
[baseMVA2, bus2, gen2, branch2, area2, gencost2] = loadcase(c);
t_is(baseMVA2,	baseMVA,	12, [t 'baseMVA']);
t_is(bus2,		bus,		12, [t 'bus']);
t_is(gen2,		gen,		12, [t 'gen']);
t_is(branch2,	branch,		12, [t 'branch']);
t_is(area2,		areas,		12, [t 'areas']);
t_is(gencost2,	gencost,	12, [t 'gencost']);

t = 'runpf(my_M_file)';
opt = mpoption('VERBOSE', 0, 'OUT_ALL', 0);
[baseMVA3, bus3, gen3, branch3, success, et] = runpf(casefile, opt);
t_ok( success, t );

t = 'runpf(my_struct)';
[baseMVA4, bus4, gen4, branch4, success, et] = runpf(c, opt);
t_ok( success, t );

t = 'runpf result comparison : ';
t_is(baseMVA3,	baseMVA4,	12, [t 'baseMVA']);
t_is(bus3,		bus4,		12, [t 'bus']);
t_is(gen3,		gen4,		12, [t 'gen']);
t_is(branch3,	branch4,	12, [t 'branch']);

t = 'runpf(modified_struct)';
c.gen(3,2) = c.gen(3,2) + 1;			%% increase gen 3 output by 1
[baseMVA5, bus5, gen5, branch5, success, et] = runpf(c, opt);
t_is(gen5(1,2), gen4(1,2) - 1, 1, t);	%% slack bus output should decrease by 1

t_end;

return;
