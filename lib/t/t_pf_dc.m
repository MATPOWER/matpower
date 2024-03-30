function t_pf_dc(quiet)
% t_pf_dc - Tests for legacy DC power flow solver.

%   MATPOWER
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

t_begin(14, quiet);

casefile = 't_case9_pf';
if quiet
    verbose = 0;
else
    verbose = 1;
end
if have_feature('octave')
    if have_feature('octave', 'vnum') >= 4
        file_in_path_warn_id = 'Octave:data-file-in-path';
    else
        file_in_path_warn_id = 'Octave:load-file-in-path';
    end
    s1 = warning('query', file_in_path_warn_id);
    warning('off', file_in_path_warn_id);
end
mpopt = mpoption('out.all', 0, 'pf.tol', 1e-9, 'verbose', 0);
mpopt = mpoption(mpopt, 'verbose', verbose);

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% network with islands
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
mpc1 = mpc;

%%-----  DC power flow  -----
%% get solved AC power flow case from MAT-file
load soln9_dcpf;        %% defines bus_soln, gen_soln, branch_soln

%% run DC PF
t = 'DC PF : ';
[baseMVA, bus, gen, branch, success, et] = rundcpf(casefile, mpopt);
t_ok(success, [t 'success']);
t_is(bus, bus_soln, 6, [t 'bus']);
t_is(gen, gen_soln, 6, [t 'gen']);
t_is(branch, branch_soln, 6, [t 'branch']);
r = rundcpf(casefile, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.bus, bus_soln, 6, [t 'bus']);
t_is(r.gen, gen_soln, 6, [t 'gen']);
t_is(r.branch, branch_soln, 6, [t 'branch']);

%% network with islands
t = sprintf('DC PF - network w/islands : ');
mpc  = mpc1;
%mpopt = mpoption(mpopt, 'out.bus', 1, 'out.gen', 1, 'out.all', -1, 'verbose', 2);
r = rundcpf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.bus( 1:9,  VA), bus_soln(:, VA), 8, [t 'voltage angles 1']);
t_is(r.bus(10:18, VA), bus_soln(:, VA), 8, [t 'voltage angles 2']);
Pg = [gen_soln(1, PG)-30; 30; gen_soln(2:3, PG)];
t_is(r.gen(1:4, PG), Pg, 8, [t 'active power generation 1']);
t_is(r.gen(5:8, PG), Pg, 8, [t 'active power generation 1']);

%% island without slack bus (catch singluar matrix?)
t = sprintf('DC PF - network w/islands w/o slack : ');
k = find(mpc.bus(:, BUS_TYPE) == REF);
mpc.bus(k(2), BUS_TYPE) = PV;
warn_state = warning;
warning('off', 'all');  %% turn of (near-)singular matrix warnings
r = rundcpf(mpc, mpopt);
warning(warn_state);
t_is(r.success, 0, 12, [t 'success = 0']);

t_end;

if have_feature('octave')
    warning(s1.state, file_in_path_warn_id);
end
