function t_opf_mips_sc(quiet)
%T_OPF_MIPS_SC  Tests for step-controlled MIPS-based AC optimal power flow.

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

num_tests = 88;

t_begin(num_tests, quiet);

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

casefile = 't_case9_opf';
if quiet
    verbose = 0;
else
    verbose = 0;
end
if have_fcn('octave')
	s1 = warning('query', 'Octave:load-file-in-path');
    warning('off', 'Octave:load-file-in-path');
end

t0 = 'MIPS-sc : ';
mpopt = mpoption('OPF_VIOLATION', 1e-6, 'PDIPM_GRADTOL', 1e-8, ...
        'PDIPM_COMPTOL', 1e-8, 'PDIPM_COSTTOL', 1e-9);
mpopt = mpoption(mpopt, 'OUT_ALL', 0, 'VERBOSE', verbose, 'OPF_ALG', 565);

%% set up indices
ib_data     = [1:BUS_AREA BASE_KV:VMIN];
ib_voltage  = [VM VA];
ib_lam      = [LAM_P LAM_Q];
ib_mu       = [MU_VMAX MU_VMIN];
ig_data     = [GEN_BUS QMAX QMIN MBASE:APF];
ig_disp     = [PG QG VG];
ig_mu       = (MU_PMAX:MU_QMIN);
ibr_data    = (1:ANGMAX);
ibr_flow    = (PF:QT);
ibr_mu      = [MU_SF MU_ST];
ibr_angmu   = [MU_ANGMIN MU_ANGMAX];

%% get solved AC power flow case from MAT-file
load soln9_opf;     %% defines bus_soln, gen_soln, branch_soln, f_soln

%% run OPF
t = t0;
[baseMVA, bus, gen, gencost, branch, f, success, et] = runopf(casefile, mpopt);
t_ok(success, [t 'success']);
t_is(f, f_soln, 3, [t 'f']);
t_is(   bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
t_is(   bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
t_is(   bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
t_is(   bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
t_is(   gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
t_is(   gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
t_is(   gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
t_is(branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
t_is(branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
t_is(branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);

%% get solved AC power flow case from MAT-file
load soln9_opf_Plim;       %% defines bus_soln, gen_soln, branch_soln, f_soln

%% run OPF with active power line limits
t = [t0 '(P line lim) : '];
mpopt1 = mpoption(mpopt, 'OPF_FLOW_LIM', 1);
[baseMVA, bus, gen, gencost, branch, f, success, et] = runopf(casefile, mpopt1);
t_ok(success, [t 'success']);
t_is(f, f_soln, 3, [t 'f']);
t_is(   bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
t_is(   bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
t_is(   bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
t_is(   bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
t_is(   gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
t_is(   gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
t_is(   gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
t_is(branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
t_is(branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
t_is(branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);

%%-----  test OPF with quadratic gen costs moved to generalized costs  -----
mpc = loadcase(casefile);
mpc.gencost = [
	2   1500    0   3   0.11    5   0;
	2   2000    0   3   0.085   1.2 0;
	2   3000    0   3   0.1225  1   0;
];
[baseMVA, bus_soln, gen_soln, gencost, branch_soln, f_soln, success, et] = runopf(mpc, mpopt);
branch_soln = branch_soln(:,1:MU_ST);

A = sparse(0,0);
l = [];
u = [];
nb = size(mpc.bus, 1);      % number of buses
ng = size(mpc.gen, 1);      % number of gens
thbas = 1;                thend    = thbas+nb-1;
vbas     = thend+1;       vend     = vbas+nb-1;
pgbas    = vend+1;        pgend    = pgbas+ng-1;
qgbas    = pgend+1;       qgend    = qgbas+ng-1;
nxyz = 2*nb + 2*ng;
N = sparse((1:ng)', (pgbas:pgend)', mpc.baseMVA * ones(ng,1), ng, nxyz);
fparm = ones(ng,1) * [ 1 0 0 1 ];
[junk, ix] = sort(mpc.gen(:, 1));
H = 2 * spdiags(mpc.gencost(ix, 5), 0, ng, ng);
Cw = mpc.gencost(ix, 6);
mpc.gencost(:, 5:7) = 0;

%% run OPF with quadratic gen costs moved to generalized costs
t = [t0 'w/quadratic generalized gen cost : '];
[r, success] = opf(mpc, A, l, u, mpopt, N, fparm, H, Cw);
[f, bus, gen, branch] = deal(r.f, r.bus, r.gen, r.branch);
t_ok(success, [t 'success']);
t_is(f, f_soln, 3, [t 'f']);
t_is(   bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
t_is(   bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
t_is(   bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
t_is(   bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
t_is(   gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
t_is(   gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
t_is(   gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
t_is(branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
t_is(branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
t_is(branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);
t_is(r.cost.usr, f, 12, [t 'user cost']);

%%-----  run OPF with extra linear user constraints & costs  -----
%% single new z variable constrained to be greater than or equal to
%% deviation from 1 pu voltage at bus 1, linear cost on this z
%% get solved AC power flow case from MAT-file
load soln9_opf_extras1;   %% defines bus_soln, gen_soln, branch_soln, f_soln
A = sparse([1;1;2;2],[10;25;10;25],[-1;1;1;1],2,25);
u = [Inf; Inf];
l = [-1; 1];

N = sparse(1, 25, 1, 1, 25);    %% new z variable only
fparm = [1 0 0 1];              %% w = r = z
H = sparse(1,1);                %% no quadratic term
Cw = 100;

t = [t0 'w/extra constraints & costs 1 : '];
[r, success] = opf(casefile, A, l, u, mpopt, N, fparm, H, Cw);
[f, bus, gen, branch] = deal(r.f, r.bus, r.gen, r.branch);
t_ok(success, [t 'success']);
t_is(f, f_soln, 3, [t 'f']);
t_is(   bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
t_is(   bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
t_is(   bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
t_is(   bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
t_is(   gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
t_is(   gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
t_is(   gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
t_is(branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
t_is(branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
t_is(branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);
t_is(r.var.val.z, 0.025419, 6, [t 'user variable']);
t_is(r.cost.usr, 2.5419, 4, [t 'user cost']);

%%-----  test OPF with capability curves  -----
mpc = loadcase('t_case9_opfv2');
%% remove angle diff limits
mpc.branch(1, ANGMAX) = 360;
mpc.branch(9, ANGMIN) = -360;

%% get solved AC power flow case from MAT-file
load soln9_opf_PQcap;   %% defines bus_soln, gen_soln, branch_soln, f_soln
	
%% run OPF with capability curves
t = [t0 'w/capability curves : '];
[baseMVA, bus, gen, gencost, branch, f, success, et] = runopf(mpc, mpopt);
t_ok(success, [t 'success']);
t_is(f, f_soln, 3, [t 'f']);
t_is(   bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
t_is(   bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
t_is(   bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
t_is(   bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
t_is(   gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
t_is(   gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
t_is(   gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
t_is(branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
t_is(branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
t_is(branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);

%%-----  test OPF with angle difference limits  -----
mpc = loadcase('t_case9_opfv2');
%% remove capability curves
mpc.gen(2:3, [PC1, PC2, QC1MIN, QC1MAX, QC2MIN, QC2MAX]) = zeros(2,6);

%% get solved AC power flow case from MAT-file
load soln9_opf_ang;   %% defines bus_soln, gen_soln, branch_soln, f_soln
	
%% run OPF with angle difference limits
t = [t0 'w/angle difference limits : '];
[baseMVA, bus, gen, gencost, branch, f, success, et] = runopf(mpc, mpopt);
t_ok(success, [t 'success']);
t_is(f, f_soln, 3, [t 'f']);
t_is(   bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
t_is(   bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
t_is(   bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
t_is(   bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  1, [t 'bus mu']);
t_is(   gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
t_is(   gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
t_is(   gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
t_is(branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
t_is(branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
t_is(branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);
t_is(branch(:,ibr_angmu ), branch_soln(:,ibr_angmu ),  2, [t 'branch angle mu']);

%%-----  test OPF with ignored angle difference limits  -----
%% get solved AC power flow case from MAT-file
load soln9_opf;   %% defines bus_soln, gen_soln, branch_soln, f_soln

%% run OPF with ignored angle difference limits
t = [t0 'w/ignored angle difference limits : '];
mpopt1 = mpoption(mpopt, 'OPF_IGNORE_ANG_LIM', 1);
[baseMVA, bus, gen, gencost, branch, f, success, et] = runopf(mpc, mpopt1);
%% ang limits are not in this solution data, so let's remove them
branch(1, ANGMAX) = 360;
branch(9, ANGMIN) = -360;
t_ok(success, [t 'success']);
t_is(f, f_soln, 3, [t 'f']);
t_is(   bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
t_is(   bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
t_is(   bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
t_is(   bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
t_is(   gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
t_is(   gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
t_is(   gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
t_is(branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
t_is(branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
t_is(branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);

if have_fcn('octave')
    warning(s1.state, 'Octave:load-file-in-path');
end

t_end;
