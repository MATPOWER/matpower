function t_load2disp(quiet)
%T_LOAD2DISP  Tests for LOAD2DISP.

%   MATPOWER
%   Copyright (c) 2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

num_tests = 72;

t_begin(num_tests, quiet);

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

casefile = 't_case9_opf';
if quiet
    verbose = 0;
else
    verbose = 0;
end
if have_fcn('octave')
    if have_fcn('octave', 'vnum') >= 4
        file_in_path_warn_id = 'Octave:data-file-in-path';
    else
        file_in_path_warn_id = 'Octave:load-file-in-path';
    end
    s1 = warning('query', file_in_path_warn_id);
    warning('off', file_in_path_warn_id);
end

mpopt = mpoption('opf.violation', 1e-6, 'mips.gradtol', 1e-8, ...
        'mips.comptol', 1e-8, 'mips.costtol', 1e-9);
mpopt = mpoption(mpopt, 'out.all', 0, 'verbose', verbose, 'opf.ac.solver', 'MIPS');

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

%% get solved AC OPF case from MAT-file
load soln9_opf;     %% defines bus_soln, gen_soln, branch_soln, f_soln

%% load2disp
mpc0 = loadcase(casefile);
mpc1 = load2disp(mpc0);
ng = size(mpc0.gen, 1);
ib = find(mpc0.bus(:, PD));
ig = ng+(1:length(ib));

t = 'mpc = load2disp(mpc) : ';
t_is(size(mpc1.gen, 1), size(mpc0.gen, 1)+length(ib), 12, [t 'number of gens']);
t_is(mpc1.bus(:, PD), 0, 12, [t 'PD is zero']);
t_is(mpc1.bus(:, QD), 0, 12, [t 'QD is zero']);
t_is(mpc1.gen(ig, PG), -mpc0.bus(ib, PD), 12, [t 'PG = -PD']);
t_is(mpc1.gen(ig, PMIN), -mpc0.bus(ib, PD), 12, [t 'PMIN = -PD']);
t_is(mpc1.gen(ig, PMAX), 0, 12, [t 'PMAX = 0']);
t_is(mpc1.gen(ig, QMIN), -mpc0.bus(ib, QD), 12, [t 'QMIN = -QD']);
t_is(mpc1.gen(ig, QMAX), 0, 12, [t 'QMAX = 0']);
t_is(mpc1.gen(ig, VG), mpc0.bus(ib, VM), 12, [t 'VG']);
t_is(mpc1.gencost(ig, MODEL), 2, 12, [t 'gencost(:, MODEL)']);
t_is(mpc1.gencost(ig, NCOST), 2, 12, [t 'gencost(:, NCOST)']);
t_is(mpc1.gencost(ig, COST), 5000, 12, [t 'gencost(:, COST)']);

t = 'mpc = load2disp(mpc, '''', idx, voll) : ';
idx = [5;7];
ig = ng+(1:length(idx));
mpc2 = load2disp(mpc0, '', idx, 1000);
t_is(size(mpc2.gen, 1), size(mpc0.gen, 1)+length(idx), 12, [t 'number of gens']);
t_is(mpc2.bus(:, PD), [zeros(8,1); 125], 12, [t 'PD']);
t_is(mpc2.bus(:, QD), [zeros(8,1);  50], 12, [t 'QD']);
t_is(mpc2.gen(ig, PG), -mpc0.bus(idx, PD), 12, [t 'PG = -PD']);
t_is(mpc2.gen(ig, PMIN), -mpc0.bus(idx, PD), 12, [t 'PMIN = -PD']);
t_is(mpc2.gen(ig, PMAX), 0, 12, [t 'PMAX = 0']);
t_is(mpc2.gen(ig, QMIN), -mpc0.bus(idx, QD), 12, [t 'QMIN = -QD']);
t_is(mpc2.gen(ig, QMAX), 0, 12, [t 'QMAX = 0']);
t_is(mpc2.gen(ig, VG), mpc0.bus(idx, VM), 12, [t 'VG']);
t_is(mpc2.gencost(ig, MODEL), 2, 12, [t 'gencost(:, MODEL)']);
t_is(mpc2.gencost(ig, NCOST), 2, 12, [t 'gencost(:, NCOST)']);
t_is(mpc2.gencost(ig, COST), 1000, 12, [t 'gencost(:, COST)']);

%% run OPF
t = 'fixed load OPF soln : ';
r0 = runopf(mpc0, mpopt);
t_ok(r0.success, [t 'success']);
t_is(r0.f, f_soln, 3, [t 'f']);
t_is(r0.bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
t_is(r0.bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
t_is(r0.bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
t_is(r0.bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
t_is(r0.gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
t_is(r0.gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
t_is(r0.gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
t_is(r0.branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
t_is(r0.branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
t_is(r0.branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);

t = 'disp load OPF soln : ';
r1 = runopf(mpc1, mpopt);
t_ok(r0.success, [t 'success']);
t_is(r0.f, f_soln, 3, [t 'f']);
t_is(r0.bus(:,ib_data   ),    bus_soln(:,ib_data   ), 10, [t 'bus data']);
t_is(r0.bus(:,ib_voltage),    bus_soln(:,ib_voltage),  3, [t 'bus voltage']);
t_is(r0.bus(:,ib_lam    ),    bus_soln(:,ib_lam    ),  3, [t 'bus lambda']);
t_is(r0.bus(:,ib_mu     ),    bus_soln(:,ib_mu     ),  2, [t 'bus mu']);
t_is(r0.gen(:,ig_data   ),    gen_soln(:,ig_data   ), 10, [t 'gen data']);
t_is(r0.gen(:,ig_disp   ),    gen_soln(:,ig_disp   ),  3, [t 'gen dispatch']);
t_is(r0.gen(:,ig_mu     ),    gen_soln(:,ig_mu     ),  3, [t 'gen mu']);
t_is(r0.branch(:,ibr_data  ), branch_soln(:,ibr_data  ), 10, [t 'branch data']);
t_is(r0.branch(:,ibr_flow  ), branch_soln(:,ibr_flow  ),  3, [t 'branch flow']);
t_is(r0.branch(:,ibr_mu    ), branch_soln(:,ibr_mu    ),  2, [t 'branch mu']);

t = 'mpc = load2disp(solved_mpc) + gentype/genfuel : ';
r0.genfuel = {'ng'; 'coal'; 'hydro'};
r0.gentype = {'GT'; 'ST'; 'HY'};
mpc3 = load2disp(r0);
ib = find(mpc0.bus(:, PD));
ig = ng+(1:length(ib));
t_is(size(mpc3.gen, 1), size(mpc0.gen, 1)+length(ib), 12, [t 'number of gens']);
t_is(mpc3.bus(:, PD), 0, 12, [t 'PD is zero']);
t_is(mpc3.bus(:, QD), 0, 12, [t 'QD is zero']);
t_is(mpc3.gen(ig, PG), -mpc0.bus(ib, PD), 12, [t 'PG = -PD']);
t_is(mpc3.gen(ig, PMIN), -mpc0.bus(ib, PD), 12, [t 'PMIN = -PD']);
t_is(mpc3.gen(ig, PMAX), 0, 12, [t 'PMAX = 0']);
t_is(mpc3.gen(ig, QMIN), -mpc0.bus(ib, QD), 12, [t 'QMIN = -QD']);
t_is(mpc3.gen(ig, QMAX), 0, 12, [t 'QMAX = 0']);
t_is(mpc3.gen(ig, VG), r0.bus(ib, VM), 12, [t 'VG']);
t_is(mpc3.gencost(ig, MODEL), 2, 12, [t 'gencost(:, MODEL)']);
t_is(mpc3.gencost(ig, NCOST), 2, 12, [t 'gencost(:, NCOST)']);
t_is(mpc3.gencost(ig, COST), 5000, 12, [t 'gencost(:, COST)']);
t_ok(isfield(mpc3, 'gentype'), [t 'results.gentype']);
t_ok(isfield(mpc3, 'genfuel'), [t 'results.genfuel']);
if isfield(mpc3, 'gentype')
    t_ok(isequal(mpc3.gentype, {r0.gentype{:}, 'DL', 'DL', 'DL'}'), [t 'results.results.gentype']);
else
    t_ok(0, [t 'results.results.gentype']);
    got = mpc3.gentype
    expected = {r0.gentype{:}, 'DL', 'DL', 'DL'}'
end
if isfield(mpc3, 'genfuel')
    t_ok(isequal(mpc3.genfuel, {r0.genfuel{:}, 'dl', 'dl', 'dl'}'), [t 'results.results.genfuel']);
else
    t_ok(0, [t 'results.results.genfuel']);
    got = mpc3.genfuel
    expected = {r0.genfuel{:}, 'dl', 'dl', 'dl'}'
end

t = 'loadshed(mpc) : ';
shed = loadshed(mpc3.gen);
t_is(size(shed), [3 1], 12, [t 'size']);
t_is(shed, [0;0;0], 12, [t 'all zeros']);

t = 'loadshed(mpc, ild) : ';
ild = [6; 5];
shed = loadshed(mpc3.gen, ild);
t_is(size(shed), [2 1], 12, [t 'size']);
t_is(shed, [0;0], 12, [t 'all zeros']);

ild = find(isload(mpc3.gen));
mpc3.gen(ild, PG) = [-80; -100; -120];

t = 'loadshed(mpc) : ';
shed = loadshed(mpc3.gen);
t_is(size(shed), [3 1], 12, [t 'size']);
t_is(shed, [10;0;5], 12, [t 'all zeros']);

t = 'loadshed(mpc, ild) : ';
ild = [6; 5];
shed = loadshed(mpc3.gen, ild);
t_is(size(shed), [2 1], 12, [t 'size']);
t_is(shed, [5;0], 12, [t 'all zeros']);

if have_fcn('octave')
    warning(s1.state, file_in_path_warn_id);
end

t_end;
