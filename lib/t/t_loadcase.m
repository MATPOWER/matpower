function t_loadcase(quiet)
%T_LOADCASE  Test that LOADCASE works with a struct as well as case file.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2004-2010 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

if nargin < 1
    quiet = 0;
end

t_begin(240, quiet);

%% compare result of loading from M-file file to result of using data matrices
casefile = 't_case9_opf';
matfile  = 't_mat9_opf';
pfcasefile = 't_case9_pf';
pfmatfile  = 't_mat9_pf';
casefilev2 = 't_case9_opfv2';
matfilev2  = 't_mat9_opfv2';
pfcasefilev2 = 't_case9_pfv2';
pfmatfilev2  = 't_mat9_pfv2';

[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% read version 1 OPF data matrices
[baseMVA, bus, gen, branch, areas, gencost] = feval(casefile);
%% save as .mat file
eval(['save ' matfile '.mat baseMVA bus gen branch areas gencost']);

%% read version 2 OPF data matrices
mpc = feval(casefilev2);
%% save as .mat file
eval(['save ' matfilev2 '.mat mpc']);

%% prepare expected matrices for v1 load
%% (missing gen cap curve & branch ang diff lims)
tmp1 = {mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.areas, mpc.gencost};
tmp2 = {mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.areas, mpc.gencost};
%% remove capability curves, angle difference limits
tmp1{3}(2:3, [PC1, PC2, QC1MIN, QC1MAX, QC2MIN, QC2MAX]) = zeros(2,6);
tmp1{4}(1, ANGMAX) = 360;
tmp1{4}(9, ANGMIN) = -360;
[baseMVA, bus, gen, branch, areas, gencost] = deal(tmp1{:});

%%-----  load OPF data into individual matrices  -----
t = 'loadcase(opf_M_file_v1) without .m extension : ';
[baseMVA1, bus1, gen1, branch1, areas1, gencost1] = loadcase(casefile);
t_is(baseMVA1,  baseMVA,    12, [t 'baseMVA']);
t_is(bus1,      bus,        12, [t 'bus']);
t_is(gen1,      gen,        12, [t 'gen']);
t_is(branch1,   branch,     12, [t 'branch']);
t_is(areas1,    areas,      12, [t 'areas']);
t_is(gencost1,  gencost,    12, [t 'gencost']);

t = 'loadcase(opf_M_file_v1) with .m extension : ';
[baseMVA1, bus1, gen1, branch1, areas1, gencost1] = loadcase([casefile '.m']);
t_is(baseMVA1,  baseMVA,    12, [t 'baseMVA']);
t_is(bus1,      bus,        12, [t 'bus']);
t_is(gen1,      gen,        12, [t 'gen']);
t_is(branch1,   branch,     12, [t 'branch']);
t_is(areas1,    areas,      12, [t 'areas']);
t_is(gencost1,  gencost,    12, [t 'gencost']);

t = 'loadcase(opf_MAT_file_v1) without .mat extension : ';
[baseMVA1, bus1, gen1, branch1, areas1, gencost1] = loadcase(matfile);
t_is(baseMVA1,  baseMVA,    12, [t 'baseMVA']);
t_is(bus1,      bus,        12, [t 'bus']);
t_is(gen1,      gen,        12, [t 'gen']);
t_is(branch1,   branch,     12, [t 'branch']);
t_is(areas1,    areas,      12, [t 'areas']);
t_is(gencost1,  gencost,    12, [t 'gencost']);

t = 'loadcase(opf_MAT_file_v1) with .mat extension : ';
[baseMVA1, bus1, gen1, branch1, areas1, gencost1] = loadcase([matfile '.mat']);
t_is(baseMVA1,  baseMVA,    12, [t 'baseMVA']);
t_is(bus1,      bus,        12, [t 'bus']);
t_is(gen1,      gen,        12, [t 'gen']);
t_is(branch1,   branch,     12, [t 'branch']);
t_is(areas1,    areas,      12, [t 'areas']);
t_is(gencost1,  gencost,    12, [t 'gencost']);

%% prepare expected matrices for v2 load
[baseMVA, bus, gen, branch, areas, gencost] = deal(tmp2{:});

t = 'loadcase(opf_M_file_v2) without .m extension : ';
[baseMVA1, bus1, gen1, branch1, areas1, gencost1] = loadcase(casefilev2);
t_is(baseMVA1,  baseMVA,    12, [t 'baseMVA']);
t_is(bus1,      bus,        12, [t 'bus']);
t_is(gen1,      gen,        12, [t 'gen']);
t_is(branch1,   branch,     12, [t 'branch']);
t_is(areas1,    areas,      12, [t 'areas']);
t_is(gencost1,  gencost,    12, [t 'gencost']);

t = 'loadcase(opf_M_file_v2) with .m extension : ';
[baseMVA1, bus1, gen1, branch1, areas1, gencost1] = loadcase([casefilev2 '.m']);
t_is(baseMVA1,  baseMVA,    12, [t 'baseMVA']);
t_is(bus1,      bus,        12, [t 'bus']);
t_is(gen1,      gen,        12, [t 'gen']);
t_is(branch1,   branch,     12, [t 'branch']);
t_is(areas1,    areas,      12, [t 'areas']);
t_is(gencost1,  gencost,    12, [t 'gencost']);

t = 'loadcase(opf_MAT_file_v2) without .mat extension : ';
[baseMVA1, bus1, gen1, branch1, areas1, gencost1] = loadcase(matfilev2);
t_is(baseMVA1,  baseMVA,    12, [t 'baseMVA']);
t_is(bus1,      bus,        12, [t 'bus']);
t_is(gen1,      gen,        12, [t 'gen']);
t_is(branch1,   branch,     12, [t 'branch']);
t_is(areas1,    areas,      12, [t 'areas']);
t_is(gencost1,  gencost,    12, [t 'gencost']);

t = 'loadcase(opf_MAT_file_v2) with .mat extension : ';
[baseMVA1, bus1, gen1, branch1, areas1, gencost1] = loadcase([matfilev2 '.mat']);
t_is(baseMVA1,  baseMVA,    12, [t 'baseMVA']);
t_is(bus1,      bus,        12, [t 'bus']);
t_is(gen1,      gen,        12, [t 'gen']);
t_is(branch1,   branch,     12, [t 'branch']);
t_is(areas1,    areas,      12, [t 'areas']);
t_is(gencost1,  gencost,    12, [t 'gencost']);

%% prepare expected matrices for v1 load
[baseMVA, bus, gen, branch, areas, gencost] = deal(tmp1{:});

t = 'loadcase(opf_struct_v1) (no version): ';
[baseMVA1, bus1, gen1, branch1, areas1, gencost1] = feval(casefile);
clear c;
c.baseMVA   = baseMVA1;
c.bus       = bus1;
c.gen       = gen1;
c.branch    = branch1;
c.areas     = areas1;
c.gencost   = gencost1;
[baseMVA2, bus2, gen2, branch2, areas2, gencost2] = loadcase(c);
t_is(baseMVA2,  baseMVA,    12, [t 'baseMVA']);
t_is(bus2,      bus,        12, [t 'bus']);
t_is(gen2,      gen,        12, [t 'gen']);
t_is(branch2,   branch,     12, [t 'branch']);
t_is(areas2,    areas,      12, [t 'areas']);
t_is(gencost2,  gencost,    12, [t 'gencost']);

t = 'loadcase(opf_struct_v1) (version=''1''): ';
c.version   = '1';
[baseMVA2, bus2, gen2, branch2, areas2, gencost2] = loadcase(c);
t_is(baseMVA2,  baseMVA,    12, [t 'baseMVA']);
t_is(bus2,      bus,        12, [t 'bus']);
t_is(gen2,      gen,        12, [t 'gen']);
t_is(branch2,   branch,     12, [t 'branch']);
t_is(areas2,    areas,      12, [t 'areas']);
t_is(gencost2,  gencost,    12, [t 'gencost']);

%% prepare expected matrices for v2 load
[baseMVA, bus, gen, branch, areas, gencost] = deal(tmp2{:});

t = 'loadcase(opf_struct_v2) (no version): ';
clear c;
c.baseMVA   = baseMVA;
c.bus       = bus;
c.gen       = gen;
c.branch    = branch;
c.areas     = areas;
c.gencost   = gencost;
[baseMVA2, bus2, gen2, branch2, areas2, gencost2] = loadcase(c);
t_is(baseMVA2,  baseMVA,    12, [t 'baseMVA']);
t_is(bus2,      bus,        12, [t 'bus']);
t_is(gen2,      gen,        12, [t 'gen']);
t_is(branch2,   branch,     12, [t 'branch']);
t_is(areas2,    areas,      12, [t 'areas']);
t_is(gencost2,  gencost,    12, [t 'gencost']);

t = 'loadcase(opf_struct_v2) (version=''2''): ';
clear c;
c.baseMVA   = baseMVA;
c.bus       = bus;
c.gen       = gen;
c.branch    = branch;
c.areas     = areas;
c.gencost   = gencost;
c.version   = '2';
[baseMVA2, bus2, gen2, branch2, areas2, gencost2] = loadcase(c);
t_is(baseMVA2,  baseMVA,    12, [t 'baseMVA']);
t_is(bus2,      bus,        12, [t 'bus']);
t_is(gen2,      gen,        12, [t 'gen']);
t_is(branch2,   branch,     12, [t 'branch']);
t_is(areas2,    areas,      12, [t 'areas']);
t_is(gencost2,  gencost,    12, [t 'gencost']);

%%-----  load OPF data into struct  -----
%% prepare expected matrices for v1 load
[baseMVA, bus, gen, branch, areas, gencost] = deal(tmp1{:});

t = 'mpc = loadcase(opf_M_file_v1) without .m extension : ';
mpc1 = loadcase(casefile);
t_is(mpc1.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc1.bus,      bus,        12, [t 'bus']);
t_is(mpc1.gen,      gen,        12, [t 'gen']);
t_is(mpc1.branch,   branch,     12, [t 'branch']);
t_is(mpc1.areas,    areas,      12, [t 'areas']);
t_is(mpc1.gencost,  gencost,    12, [t 'gencost']);

t = 'mpc = loadcase(opf_M_file_v1) with .m extension : ';
mpc1 = loadcase([casefile '.m']);
t_is(mpc1.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc1.bus,      bus,        12, [t 'bus']);
t_is(mpc1.gen,      gen,        12, [t 'gen']);
t_is(mpc1.branch,   branch,     12, [t 'branch']);
t_is(mpc1.areas,    areas,      12, [t 'areas']);
t_is(mpc1.gencost,  gencost,    12, [t 'gencost']);

t = 'mpc = loadcase(opf_MAT_file_v1) without .mat extension : ';
mpc1 = loadcase(matfile);
t_is(mpc1.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc1.bus,      bus,        12, [t 'bus']);
t_is(mpc1.gen,      gen,        12, [t 'gen']);
t_is(mpc1.branch,   branch,     12, [t 'branch']);
t_is(mpc1.areas,    areas,      12, [t 'areas']);
t_is(mpc1.gencost,  gencost,    12, [t 'gencost']);

t = 'mpc = loadcase(opf_MAT_file_v1) with .mat extension : ';
mpc1 = loadcase([matfile '.mat']);
t_is(mpc1.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc1.bus,      bus,        12, [t 'bus']);
t_is(mpc1.gen,      gen,        12, [t 'gen']);
t_is(mpc1.branch,   branch,     12, [t 'branch']);
t_is(mpc1.areas,    areas,      12, [t 'areas']);
t_is(mpc1.gencost,  gencost,    12, [t 'gencost']);

%% prepare expected matrices for v2 load
[baseMVA, bus, gen, branch, areas, gencost] = deal(tmp2{:});

t = 'mpc = loadcase(opf_M_file_v2) without .m extension : ';
mpc1 = loadcase(casefilev2);
t_is(mpc1.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc1.bus,      bus,        12, [t 'bus']);
t_is(mpc1.gen,      gen,        12, [t 'gen']);
t_is(mpc1.branch,   branch,     12, [t 'branch']);
t_is(mpc1.areas,    areas,      12, [t 'areas']);
t_is(mpc1.gencost,  gencost,    12, [t 'gencost']);

t = 'mpc = loadcase(opf_M_file_v2) with .m extension : ';
mpc1 = loadcase([casefilev2 '.m']);
t_is(mpc1.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc1.bus,      bus,        12, [t 'bus']);
t_is(mpc1.gen,      gen,        12, [t 'gen']);
t_is(mpc1.branch,   branch,     12, [t 'branch']);
t_is(mpc1.areas,    areas,      12, [t 'areas']);
t_is(mpc1.gencost,  gencost,    12, [t 'gencost']);

t = 'mpc = loadcase(opf_MAT_file_v2) without .mat extension : ';
mpc1 = loadcase(matfilev2);
t_is(mpc1.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc1.bus,      bus,        12, [t 'bus']);
t_is(mpc1.gen,      gen,        12, [t 'gen']);
t_is(mpc1.branch,   branch,     12, [t 'branch']);
t_is(mpc1.areas,    areas,      12, [t 'areas']);
t_is(mpc1.gencost,  gencost,    12, [t 'gencost']);

t = 'mpc = loadcase(opf_MAT_file_v2) with .mat extension : ';
mpc1 = loadcase([matfilev2 '.mat']);
t_is(mpc1.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc1.bus,      bus,        12, [t 'bus']);
t_is(mpc1.gen,      gen,        12, [t 'gen']);
t_is(mpc1.branch,   branch,     12, [t 'branch']);
t_is(mpc1.areas,    areas,      12, [t 'areas']);
t_is(mpc1.gencost,  gencost,    12, [t 'gencost']);

%% prepare expected matrices for v1 load
[baseMVA, bus, gen, branch, areas, gencost] = deal(tmp1{:});

t = 'mpc = loadcase(opf_struct_v1) (no version): ';
[baseMVA1, bus1, gen1, branch1, areas1, gencost1] = feval(casefile);
clear c;
c.baseMVA   = baseMVA1;
c.bus       = bus1;
c.gen       = gen1;
c.branch    = branch1;
c.areas     = areas1;
c.gencost   = gencost1;
mpc2 = loadcase(c);
t_is(mpc2.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc2.bus,      bus,        12, [t 'bus']);
t_is(mpc2.gen,      gen,        12, [t 'gen']);
t_is(mpc2.branch,   branch,     12, [t 'branch']);
t_is(mpc2.areas,    areas,      12, [t 'areas']);
t_is(mpc2.gencost,  gencost,    12, [t 'gencost']);

t = 'mpc = loadcase(opf_struct_v1) (version=''1''): ';
c.version   = '1';
mpc2 = loadcase(c);
t_is(mpc2.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc2.bus,      bus,        12, [t 'bus']);
t_is(mpc2.gen,      gen,        12, [t 'gen']);
t_is(mpc2.branch,   branch,     12, [t 'branch']);
t_is(mpc2.areas,    areas,      12, [t 'areas']);
t_is(mpc2.gencost,  gencost,    12, [t 'gencost']);

%% prepare expected matrices for v2 load
[baseMVA, bus, gen, branch, areas, gencost] = deal(tmp2{:});

t = 'mpc = loadcase(opf_struct_v2) (no version): ';
clear c;
c.baseMVA   = baseMVA;
c.bus       = bus;
c.gen       = gen;
c.branch    = branch;
c.areas     = areas;
c.gencost   = gencost;
mpc2 = loadcase(c);
t_is(mpc2.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc2.bus,      bus,        12, [t 'bus']);
t_is(mpc2.gen,      gen,        12, [t 'gen']);
t_is(mpc2.branch,   branch,     12, [t 'branch']);
t_is(mpc2.areas,    areas,      12, [t 'areas']);
t_is(mpc2.gencost,  gencost,    12, [t 'gencost']);
t_ok(strcmp(mpc2.version, '2'), [t 'version']);  

t = 'mpc = loadcase(opf_struct_v2) (version=''2''): ';
clear c;
c.baseMVA   = baseMVA;
c.bus       = bus;
c.gen       = gen;
c.branch    = branch;
c.areas     = areas;
c.gencost   = gencost;
c.version   = '2';
mpc2 = loadcase(c);
t_is(mpc2.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc2.bus,      bus,        12, [t 'bus']);
t_is(mpc2.gen,      gen,        12, [t 'gen']);
t_is(mpc2.branch,   branch,     12, [t 'branch']);
t_is(mpc2.areas,    areas,      12, [t 'areas']);
t_is(mpc2.gencost,  gencost,    12, [t 'gencost']);


%% read version 1 PF data matrices
[baseMVA, bus, gen, branch] = feval(pfcasefile);
eval(['save ' pfmatfile '.mat baseMVA bus gen branch']);

%% read version 2 PF data matrices
mpc = feval(pfcasefilev2);
tmp = {mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch};
[baseMVA, bus, gen, branch] = deal(tmp{:});
%% save as .mat file
eval(['save ' pfmatfilev2 '.mat mpc']);

%%-----  load PF data into individual matrices  -----
t = 'loadcase(pf_M_file_v1) without .m extension : ';
[baseMVA1, bus1, gen1, branch1] = loadcase(pfcasefile);
t_is(baseMVA1,  baseMVA,    12, [t 'baseMVA']);
t_is(bus1,      bus,        12, [t 'bus']);
t_is(gen1,      gen,        12, [t 'gen']);
t_is(branch1,   branch,     12, [t 'branch']);

t = 'loadcase(pf_M_file_v1) with .m extension : ';
[baseMVA1, bus1, gen1, branch1] = loadcase([pfcasefile '.m']);
t_is(baseMVA1,  baseMVA,    12, [t 'baseMVA']);
t_is(bus1,      bus,        12, [t 'bus']);
t_is(gen1,      gen,        12, [t 'gen']);
t_is(branch1,   branch,     12, [t 'branch']);

t = 'loadcase(pf_MAT_file_v1) without .mat extension : ';
[baseMVA1, bus1, gen1, branch1] = loadcase(pfmatfile);
t_is(baseMVA1,  baseMVA,    12, [t 'baseMVA']);
t_is(bus1,      bus,        12, [t 'bus']);
t_is(gen1,      gen,        12, [t 'gen']);
t_is(branch1,   branch,     12, [t 'branch']);

t = 'loadcase(pf_MAT_file_v1) with .mat extension : ';
[baseMVA1, bus1, gen1, branch1] = loadcase([pfmatfile '.mat']);
t_is(baseMVA1,  baseMVA,    12, [t 'baseMVA']);
t_is(bus1,      bus,        12, [t 'bus']);
t_is(gen1,      gen,        12, [t 'gen']);
t_is(branch1,   branch,     12, [t 'branch']);

t = 'loadcase(pf_M_file_v2) without .m extension : ';
[baseMVA1, bus1, gen1, branch1] = loadcase(pfcasefilev2);
t_is(baseMVA1,  baseMVA,    12, [t 'baseMVA']);
t_is(bus1,      bus,        12, [t 'bus']);
t_is(gen1,      gen,        12, [t 'gen']);
t_is(branch1,   branch,     12, [t 'branch']);

t = 'loadcase(pf_M_file_v2) with .m extension : ';
[baseMVA1, bus1, gen1, branch1] = loadcase([pfcasefilev2 '.m']);
t_is(baseMVA1,  baseMVA,    12, [t 'baseMVA']);
t_is(bus1,      bus,        12, [t 'bus']);
t_is(gen1,      gen,        12, [t 'gen']);
t_is(branch1,   branch,     12, [t 'branch']);

t = 'loadcase(pf_MAT_file_v2) without .mat extension : ';
[baseMVA1, bus1, gen1, branch1] = loadcase(pfmatfilev2);
t_is(baseMVA1,  baseMVA,    12, [t 'baseMVA']);
t_is(bus1,      bus,        12, [t 'bus']);
t_is(gen1,      gen,        12, [t 'gen']);
t_is(branch1,   branch,     12, [t 'branch']);

t = 'loadcase(pf_MAT_file_v2) with .mat extension : ';
[baseMVA1, bus1, gen1, branch1] = loadcase([pfmatfilev2 '.mat']);
t_is(baseMVA1,  baseMVA,    12, [t 'baseMVA']);
t_is(bus1,      bus,        12, [t 'bus']);
t_is(gen1,      gen,        12, [t 'gen']);
t_is(branch1,   branch,     12, [t 'branch']);

t = 'loadcase(pf_struct_v1) (no version): ';
[baseMVA1, bus1, gen1, branch1] = feval(pfcasefile);
clear c;
c.baseMVA   = baseMVA1;
c.bus       = bus1;
c.gen       = gen1;
c.branch    = branch1;
[baseMVA2, bus2, gen2, branch2] = loadcase(c);
t_is(baseMVA2,  baseMVA,    12, [t 'baseMVA']);
t_is(bus2,      bus,        12, [t 'bus']);
t_is(gen2,      gen,        12, [t 'gen']);
t_is(branch2,   branch,     12, [t 'branch']);

t = 'loadcase(pf_struct_v1) (version=''1''): ';
c.version   = '1';
[baseMVA2, bus2, gen2, branch2] = loadcase(c);
t_is(baseMVA2,  baseMVA,    12, [t 'baseMVA']);
t_is(bus2,      bus,        12, [t 'bus']);
t_is(gen2,      gen,        12, [t 'gen']);
t_is(branch2,   branch,     12, [t 'branch']);

t = 'loadcase(pf_struct_v2) : ';
clear c;
c.baseMVA   = baseMVA;
c.bus       = bus;
c.gen       = gen;
c.branch    = branch;
c.version   = '2';
[baseMVA2, bus2, gen2, branch2] = loadcase(c);
t_is(baseMVA2,  baseMVA,    12, [t 'baseMVA']);
t_is(bus2,      bus,        12, [t 'bus']);
t_is(gen2,      gen,        12, [t 'gen']);
t_is(branch2,   branch,     12, [t 'branch']);












%%-----  load PF data into struct  -----
t = 'mpc = loadcase(pf_M_file_v1) without .m extension : ';
mpc1 = loadcase(pfcasefile);
t_is(mpc1.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc1.bus,      bus,        12, [t 'bus']);
t_is(mpc1.gen,      gen,        12, [t 'gen']);
t_is(mpc1.branch,   branch,     12, [t 'branch']);

t = 'mpc = loadcase(pf_M_file_v1) with .m extension : ';
mpc1 = loadcase([pfcasefile '.m']);
t_is(mpc1.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc1.bus,      bus,        12, [t 'bus']);
t_is(mpc1.gen,      gen,        12, [t 'gen']);
t_is(mpc1.branch,   branch,     12, [t 'branch']);

t = 'mpc = loadcase(pf_MAT_file_v1) without .mat extension : ';
mpc1 = loadcase(pfmatfile);
t_is(mpc1.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc1.bus,      bus,        12, [t 'bus']);
t_is(mpc1.gen,      gen,        12, [t 'gen']);
t_is(mpc1.branch,   branch,     12, [t 'branch']);

t = 'mpc = loadcase(pf_MAT_file_v1) with .mat extension : ';
mpc1 = loadcase([pfmatfile '.mat']);
t_is(mpc1.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc1.bus,      bus,        12, [t 'bus']);
t_is(mpc1.gen,      gen,        12, [t 'gen']);
t_is(mpc1.branch,   branch,     12, [t 'branch']);

t = 'mpc = loadcase(pf_M_file_v2) without .m extension : ';
mpc1 = loadcase(pfcasefilev2);
t_is(mpc1.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc1.bus,      bus,        12, [t 'bus']);
t_is(mpc1.gen,      gen,        12, [t 'gen']);
t_is(mpc1.branch,   branch,     12, [t 'branch']);

t = 'mpc = loadcase(pf_M_file_v2) with .m extension : ';
mpc1 = loadcase([pfcasefilev2 '.m']);
t_is(mpc1.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc1.bus,      bus,        12, [t 'bus']);
t_is(mpc1.gen,      gen,        12, [t 'gen']);
t_is(mpc1.branch,   branch,     12, [t 'branch']);

t = 'mpc = loadcase(pf_MAT_file_v2) without .mat extension : ';
mpc1 = loadcase(pfmatfilev2);
t_is(mpc1.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc1.bus,      bus,        12, [t 'bus']);
t_is(mpc1.gen,      gen,        12, [t 'gen']);
t_is(mpc1.branch,   branch,     12, [t 'branch']);

t = 'mpc = loadcase(pf_MAT_file_v2) with .mat extension : ';
mpc1 = loadcase([pfmatfilev2 '.mat']);
t_is(mpc1.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc1.bus,      bus,        12, [t 'bus']);
t_is(mpc1.gen,      gen,        12, [t 'gen']);
t_is(mpc1.branch,   branch,     12, [t 'branch']);

t = 'mpc = loadcase(pf_struct_v1) (no version): ';
[baseMVA1, bus1, gen1, branch1] = feval(pfcasefile);
clear c;
c.baseMVA   = baseMVA1;
c.bus       = bus1;
c.gen       = gen1;
c.branch    = branch1;
mpc2 = loadcase(c);
t_is(mpc2.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc2.bus,      bus,        12, [t 'bus']);
t_is(mpc2.gen,      gen,        12, [t 'gen']);
t_is(mpc2.branch,   branch,     12, [t 'branch']);

t = 'mpc = loadcase(pf_struct_v1) (version=''1''): ';
c.version   = '1';
mpc2 = loadcase(c);
t_is(mpc2.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc2.bus,      bus,        12, [t 'bus']);
t_is(mpc2.gen,      gen,        12, [t 'gen']);
t_is(mpc2.branch,   branch,     12, [t 'branch']);

t = 'mpc = loadcase(pf_struct_v2) : ';
clear c;
c.baseMVA   = baseMVA;
c.bus       = bus;
c.gen       = gen;
c.branch    = branch;
c.version   = '2';
mpc2 = loadcase(c);
t_is(mpc2.baseMVA,  baseMVA,    12, [t 'baseMVA']);
t_is(mpc2.bus,      bus,        12, [t 'bus']);
t_is(mpc2.gen,      gen,        12, [t 'gen']);
t_is(mpc2.branch,   branch,     12, [t 'branch']);

%% cleanup
delete([ matfile '.mat' ]);
delete([ pfmatfile '.mat' ]);
delete([ matfilev2 '.mat' ]);
delete([ pfmatfilev2 '.mat' ]);

t = 'runpf(my_M_file)';
mpopt = mpoption('verbose', 0, 'out.all', 0);
[baseMVA3, bus3, gen3, branch3, success, et] = runpf(pfcasefile, mpopt);
t_ok( success, t );

t = 'runpf(my_struct)';
[baseMVA4, bus4, gen4, branch4, success, et] = runpf(c, mpopt);
t_ok( success, t );

t = 'runpf result comparison : ';
t_is(baseMVA3,  baseMVA4,   12, [t 'baseMVA']);
t_is(bus3,      bus4,       12, [t 'bus']);
t_is(gen3,      gen4,       12, [t 'gen']);
t_is(branch3,   branch4,    12, [t 'branch']);

t = 'runpf(modified_struct)';
c.gen(3,2) = c.gen(3,2) + 1;            %% increase gen 3 output by 1
[baseMVA5, bus5, gen5, branch5, success, et] = runpf(c, mpopt);
t_is(gen5(1,2), gen4(1,2) - 1, 1, t);   %% slack bus output should decrease by 1

t_end;
