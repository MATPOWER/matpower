function t_pf(quiet)
%T_PF  Tests for power flow solvers.

%   MATPOWER
%   Copyright (c) 2004-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

t_begin(62, quiet);

casefile = 't_case9_pf';
if quiet
    verbose = 0;
else
    verbose = 1;
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
mpopt = mpoption('out.all', 0, 'verbose', verbose);

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% get solved AC power flow case from MAT-file
load soln9_pf;      %% defines bus_soln, gen_soln, branch_soln

%% run Newton PF
t = 'Newton PF : ';
mpopt = mpoption(mpopt, 'pf.alg', 'NR');
[baseMVA, bus, gen, branch, success, et] = runpf(casefile, mpopt);
t_ok(success, [t 'success']);
t_is(bus, bus_soln, 6, [t 'bus']);
t_is(gen, gen_soln, 6, [t 'gen']);
t_is(branch, branch_soln, 6, [t 'branch']);

%% run fast-decoupled PF (XB version)
t = 'Fast Decoupled (XB) PF : ';
mpopt = mpoption(mpopt, 'pf.alg', 'FDXB');
[baseMVA, bus, gen, branch, success, et] = runpf(casefile, mpopt);
t_ok(success, [t 'success']);
t_is(bus, bus_soln, 6, [t 'bus']);
t_is(gen, gen_soln, 6, [t 'gen']);
t_is(branch, branch_soln, 6, [t 'branch']);

%% run fast-decoupled PF (BX version)
t = 'Fast Decoupled (BX) PF : ';
mpopt = mpoption(mpopt, 'pf.alg', 'FDBX');
[baseMVA, bus, gen, branch, success, et] = runpf(casefile, mpopt);
t_ok(success, [t 'success']);
t_is(bus, bus_soln, 6, [t 'bus']);
t_is(gen, gen_soln, 6, [t 'gen']);
t_is(branch, branch_soln, 6, [t 'branch']);

%% run Gauss-Seidel PF
t = 'Gauss-Seidel PF : ';
mpopt = mpoption(mpopt, 'pf.alg', 'GS');
[baseMVA, bus, gen, branch, success, et] = runpf(casefile, mpopt);
t_ok(success, [t 'success']);
t_is(bus, bus_soln, 5, [t 'bus']);
t_is(gen, gen_soln, 5, [t 'gen']);
t_is(branch, branch_soln, 5, [t 'branch']);

%% get solved AC power flow case from MAT-file
load soln9_dcpf;        %% defines bus_soln, gen_soln, branch_soln

%% run DC PF
t = 'DC PF : ';
[baseMVA, bus, gen, branch, success, et] = rundcpf(casefile, mpopt);
t_ok(success, [t 'success']);
t_is(bus, bus_soln, 6, [t 'bus']);
t_is(gen, gen_soln, 6, [t 'gen']);
t_is(branch, branch_soln, 6, [t 'branch']);

%% check Qg distribution, when Qmin = Qmax
t = 'check Qg : ';
mpopt = mpoption(mpopt, 'pf.alg', 'NR', 'verbose', 0);
mpc = loadcase(casefile);
mpc.gen(1, [QMIN QMAX]) = [20 20];
[baseMVA, bus, gen, branch, success, et] = runpf(mpc, mpopt);
t_ok(success, [t 'success']);
t_is(gen(1, QG), 24.07, 2, [t 'single gen, Qmin = Qmax']);

mpc.gen = [mpc.gen(1, :); mpc.gen];
mpc.gen(1, [QMIN QMAX]) = [10 10];
mpc.gen(2, [QMIN QMAX]) = [0 50];
[baseMVA, bus, gen, branch, success, et] = runpf(mpc, mpopt);
t_ok(success, [t 'success']);
t_is(gen(1:2, QG), [10; 14.07], 2, [t '2 gens, Qmin = Qmax for one']);

mpc.gen(1, [QMIN QMAX]) = [10 10];
mpc.gen(2, [QMIN QMAX]) = [-50 -50];
[baseMVA, bus, gen, branch, success, et] = runpf(mpc, mpopt);
t_ok(success, [t 'success']);
t_is(gen(1:2, QG), [42.03; -17.97], 2, [t '2 gens, Qmin = Qmax for both']);

mpc.gen(1, [QMIN QMAX]) = [0 50];
mpc.gen(2, [QMIN QMAX]) = [0 100];
[baseMVA, bus, gen, branch, success, et] = runpf(mpc, mpopt);
t_ok(success, [t 'success']);
t_is(gen(1:2, QG), [8.02; 16.05], 2, [t '2 gens, proportional']);

mpc.gen(1, [QMIN QMAX]) = [-50 0];
mpc.gen(2, [QMIN QMAX]) = [50 150];
[baseMVA, bus, gen, branch, success, et] = runpf(mpc, mpopt);
t_ok(success, [t 'success']);
t_is(gen(1:2, QG), [-50+8.02; 50+16.05], 2, [t '2 gens, proportional']);

mpc.gen(1, [QMIN QMAX]) = [-50 Inf];
mpc.gen(2, [QMIN QMAX]) = [50 150];
[baseMVA, bus, gen, branch, success, et] = runpf(mpc, mpopt);
t_ok(success, [t 'success']);
t_is(gen(1:2, QG), [-31.61; 55.68], 2, [t '2 gens, one infinite range']);

mpc.gen(1, [QMIN QMAX]) = [-50 Inf];
mpc.gen(2, [QMIN QMAX]) = [50 Inf];
[baseMVA, bus, gen, branch, success, et] = runpf(mpc, mpopt);
t_ok(success, [t 'success']);
t_is(gen(1:2, QG), [-33.12; 57.18], 2, [t '2 gens, both infinite range']);

mpc.gen(1, [QMIN QMAX]) = [-50 Inf];
mpc.gen(2, [QMIN QMAX]) = [-Inf 150];
[baseMVA, bus, gen, branch, success, et] = runpf(mpc, mpopt);
t_ok(success, [t 'success']);
t_is(gen(1:2, QG), [76.07; -52], 2, [t '2 gens, both infinite range']);

mpc.gen(1, [QMIN QMAX]) = [-Inf Inf];
mpc.gen(2, [QMIN QMAX]) = [-Inf Inf];
[baseMVA, bus, gen, branch, success, et] = runpf(mpc, mpopt);
t_ok(success, [t 'success']);
t_is(gen(1:2, QG), [12.03; 12.03], 2, [t '2 gens, both infinite range']);

t = 'reactive generation allocation : ';
mpc = loadcase(casefile);
%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	0	0	300	-300	1	100	1	250	10	0	0	0	0	0	0	0	0	0	0	0;
	2	54	0	0	-5	1	100	1	300	10	0	0	0	0	0	0	0	0	0	0	0;
	2	54	0	5	-5	1	100	1	300	10	0	0	0	0	0	0	0	0	0	0	0;
	2	55	0	25	 10	1	100	1	300	10	0	0	0	0	0	0	0	0	0	0	0;
	30	25	1	300	-300	1	100	1	270	10	0	0	0	0	0	0	0	0	0	0	0;
	30	30	2	300	-300	1	100	1	270	10	0	0	0	0	0	0	0	0	0	0	0;
	30	30	-3	300	-300	1	100	1	270	10	0	0	0	0	0	0	0	0	0	0	0;
];
mpc.bus(3, BUS_TYPE) = PQ;
r = runpf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.gen(2:4, QG), [-5; -5; 10] + [1; 2; 3]*1.989129794, 8, [t 'PV bus']);
t_is(r.gen(5:7, QG), [1; 2; -3], 8, [t 'PQ bus']);

%% network with islands
t = 'network w/islands : DC PF : ';
mpc0 = loadcase(casefile);
mpc0.gen(1, PG) = 60;
mpc0.gen(1, [PMIN PMAX QMIN QMAX PG QG]) = mpc0.gen(1, [PMIN PMAX QMIN QMAX PG QG]) / 2;
mpc0.gen = [mpc0.gen(1, :); mpc0.gen];
mpc1 = mpc0;
mpc  = mpc0;
nb = size(mpc.bus, 1);
mpc1.bus(:, BUS_I)      = mpc1.bus(:, BUS_I) + nb;
mpc1.branch(:, F_BUS)   = mpc1.branch(:, F_BUS) + nb;
mpc1.branch(:, T_BUS)   = mpc1.branch(:, T_BUS) + nb;
mpc1.gen(:, GEN_BUS)    = mpc1.gen(:, GEN_BUS) + nb;
mpc.bus         = [mpc.bus; mpc1.bus];
mpc.branch      = [mpc.branch; mpc1.branch];
mpc.gen         = [mpc.gen; mpc1.gen];
%mpopt = mpoption(mpopt, 'out.bus', 1, 'out.gen', 1, 'out.all', -1, 'verbose', 2);
mpopt = mpoption(mpopt, 'verbose', verbose);
r = rundcpf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.bus( 1:9,  VA), bus_soln(:, VA), 8, [t 'voltage angles 1']);
t_is(r.bus(10:18, VA), bus_soln(:, VA), 8, [t 'voltage angles 2']);
Pg = [gen_soln(1, PG)-30; 30; gen_soln(2:3, PG)];
t_is(r.gen(1:4, PG), Pg, 8, [t 'active power generation 1']);
t_is(r.gen(5:8, PG), Pg, 8, [t 'active power generation 1']);

t = 'network w/islands : AC PF : ';
%% get solved AC power flow case from MAT-file
load soln9_pf;      %% defines bus_soln, gen_soln, branch_soln
r = runpf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.bus( 1:9,  VA), bus_soln(:, VA), 8, [t 'voltage angles 1']);
t_is(r.bus(10:18, VA), bus_soln(:, VA), 8, [t 'voltage angles 2']);
Pg = [gen_soln(1, PG)-30; 30; gen_soln(2:3, PG)];
t_is(r.gen(1:4, PG), Pg, 8, [t 'active power generation 1']);
t_is(r.gen(5:8, PG), Pg, 8, [t 'active power generation 1']);

%% island without slack bus (catch singluar matrix?)
t = 'network w/islands w/o slack : DC PF : ';
k = find(mpc.bus(:, BUS_TYPE) == REF);
mpc.bus(k(2), BUS_TYPE) = PV;
warn_state = warning;
warning('off', 'all');  %% turn of (near-)singular matrix warnings
r = rundcpf(mpc, mpopt);
warning(warn_state);
t_ok(~r.success, [t 'success']);

t = 'all buses isolated : ';
mpc.bus(:, BUS_TYPE) = NONE;
try
    r = runpf(mpc, mpopt);
    t_is(r.success, 0, 12, [t 'success = 0']);
catch
    t_ok(0, [t 'unexpected fatal error']);
end

%% case 14 with Q limits
t = 'pf.enforce_q_lims == 0 : ';
mpc = loadcase('case14');
mpc.gen(1, QMIN) = -10;
mpc.gen(:, QMAX) = [10; 30; 29; 15; 15];
bt0 = mpc.bus(:, BUS_TYPE);
bt = bt0;
mpopt = mpoption(mpopt, 'pf.enforce_q_lims', 0);
r = runpf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.bus(:, BUS_TYPE), bt, 12, [t 'bus type']);
t_is(r.gen(:, QG), [-16.549300542; 43.557100134; 25.075348495; 12.730944405; 17.623451366], 8, [t 'Qg']);

t = 'pf.enforce_q_lims == 1 : ';
mpopt = mpoption(mpopt, 'pf.enforce_q_lims', 1);
r = runpf(mpc, mpopt);
bt = bt0;
bt([1 2 3 8]) = [PQ PQ REF PQ];
t_ok(~r.success, [t '~success']);
t_is(r.bus(:, BUS_TYPE), bt, 12, [t 'bus type']);
t_is(r.gen(:, QG), [-10; 30; 31.608422873; 16.420423190; 15], 4, [t 'Qg']);

t = 'pf.enforce_q_lims == 2 : ';
mpopt = mpoption(mpopt, 'pf.enforce_q_lims', 2);
r = runpf(mpc, mpopt);
bt = bt0;
bt([1 2 3 6 8]) = [REF PQ PQ PQ PQ];
t_ok(r.success, [t 'success']);
t_is(r.bus(:, BUS_TYPE), bt, 12, [t 'bus type']);
t_is(r.gen(:, QG), [-6.30936644; 30; 29; 15; 15], 8, [t 'Qg']);

t_end;

if have_fcn('octave')
    warning(s1.state, file_in_path_warn_id);
end
