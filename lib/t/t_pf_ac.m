function t_pf_ac(quiet)
% t_pf_ac - Tests for legacy AC power flow solvers.

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

%%  alg         name                            check       opts
cfg = {
    {'NR',      'Newton (default, power-polar)',[]          []  },
    {'NR-SP',   'Newton (power-polar)',         [],         []  },
    {'NR-SC',   'Newton (power-cartesian)',     [],         []  },
    {'NR-SH',   'Newton (power-hybrid)',        [],         []  },
    {'NR-IP',   'Newton (current-polar)',       [],         []  },
    {'NR-IC',   'Newton (current-cartesian)',   [],         []  },
    {'NR-IH',   'Newton (current-hybrid)',      [],         []  },
    {'FDXB',    'Fast Decoupled (XB)',          [],         []  },
    {'FDBX',    'Fast Decoupled (BX)',          [],         []  },
    {'GS',      'Gauss-Seidel',                 [],         []  },
    {'ZG',      'Implicit Z-bus Gauss',         [],         []  },
};
if have_feature('mp_core')
    cfg{end+1} = {'FSOLVE',  'fsolve (power-polar)',         'fsolve',   []  };
    cfg{end+1} = {'FSOLVE',  'fsolve (power-cartesian)',     'fsolve',   {'pf.v_cartesian', 1}  };
    cfg{end+1} = {'FSOLVE',  'fsolve (current-polar)',       'fsolve',   {'pf.current_balance', 1, 'pf.tol', 1e-10}  };
    cfg{end+1} = {'FSOLVE',  'fsolve (current-cartesian)',   'fsolve',   {'pf.v_cartesian', 1, 'pf.current_balance', 1, 'pf.tol', 1e-10}  };
%     cfg{end+1} = {'FSOLVE',  'fsolve (power-polar)',         'fsolve',   struct('Algorithm', 'trust-region-dogleg')  },
%     cfg{end+1} = {'FSOLVE',  'fsolve (power-polar)',         'fsolve',   struct('Algorithm', 'trust-region-reflective')  },
%     cfg{end+1} = {'FSOLVE',  'fsolve (power-polar)',         'fsolve',   struct('Algorithm', 'levenberg-marquardt', 'TolX', 1e-11) },
end
% %%  alg         name                            check       opts
% cfg = {
%     {'NR',      'Newton (default, power-polar)',[]          []  },
% };

ntests = 48;
t_begin(length(cfg)*ntests + 6, quiet);

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
mpopt0 = mpoption('out.all', 0, 'pf.tol', 1e-9, 'verbose', 0);
mpopt0 = mpoption(mpopt0, 'verbose', verbose);

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
mpc0.gen(1, QG) = 10;
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

%%-----  AC power flow  -----
%% get solved AC power flow case from MAT-file
load soln9_pf;      %% defines bus_soln, gen_soln, branch_soln
soln_vg = load('soln9_pf_vg');

%% run AC PF
for k = 1:length(cfg)
  if ~isempty(cfg{k}{3}) && ~have_fcn(cfg{k}{3})
    t_skip(ntests, sprintf('%s not available', cfg{k}{3}));
  else
    t = sprintf('AC PF - %s : ', cfg{k}{2});
    mpopt = mpoption(mpopt0, 'pf.alg', cfg{k}{1});
    if ~isempty(cfg{k}{4}) && iscell(cfg{k}{4})
        mpopt = mpoption(mpopt, cfg{k}{4}{:});
    end
    [baseMVA, bus, gen, branch, success, et] = runpf(casefile, mpopt);
    t_ok(success, [t 'success']);
    t_is(bus, bus_soln, 6, [t 'bus']);
    t_is(gen, gen_soln, 6, [t 'gen']);
    t_is(branch, branch_soln, 6, [t 'branch']);

    r = runpf(casefile, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.bus, bus_soln, 6, [t 'bus']);
    t_is(r.gen, gen_soln, 6, [t 'gen']);
    t_is(r.branch, branch_soln, 6, [t 'branch']);

    %% check when Vg ~= 1
    t = sprintf('%s - Vg ~= 1 : ', cfg{k}{1});
    mpopt = mpoption(mpopt, 'verbose', 0);
    mpc = loadcase(casefile);
    mpc.gen(:, VG) = [1.04; 1.025; 1.025];
    r = runpf(mpc, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.bus, soln_vg.bus_soln, 6, [t 'bus']);
    t_is(r.gen, soln_vg.gen_soln, 6, [t 'gen']);
    t_is(r.branch, soln_vg.branch_soln, 6, [t 'branch']);

    %% check Qg distribution, when Qmin = Qmax
    t = sprintf('%s - check Qg : ', cfg{k}{1});
    mpc = loadcase(casefile);
    mpc.gen(1, [QG QMIN QMAX]) = [10 20 20];
    r = runpf(mpc, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.gen(1, QG), 24.07, 2, [t 'single gen, Qmin = Qmax']);

    mpc.gen = [mpc.gen(1, :); mpc.gen];
    mpc.gen(1, [QMIN QMAX]) = [10 10];
    mpc.gen(2, [QMIN QMAX]) = [0 50];
    r = runpf(mpc, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.gen(1:2, QG), [10; 14.07], 2, [t '2 gens, Qmin = Qmax for one']);

    mpc.gen(1, [QMIN QMAX]) = [10 10];
    mpc.gen(2, [QMIN QMAX]) = [-50 -50];
    r = runpf(mpc, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.gen(1:2, QG), [42.03; -17.97], 2, [t '2 gens, Qmin = Qmax for both']);

    mpc.gen(1, [QMIN QMAX]) = [0 50];
    mpc.gen(2, [QMIN QMAX]) = [0 100];
    r = runpf(mpc, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.gen(1:2, QG), [8.02; 16.05], 2, [t '2 gens, proportional']);

    mpc.gen(1, [QMIN QMAX]) = [-50 0];
    mpc.gen(2, [QMIN QMAX]) = [50 150];
    r = runpf(mpc, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.gen(1:2, QG), [-50+8.02; 50+16.05], 2, [t '2 gens, proportional']);

    mpc.gen(1, [QMIN QMAX]) = [-50 Inf];
    mpc.gen(2, [QMIN QMAX]) = [50 150];
    r = runpf(mpc, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.gen(1:2, QG), [-31.61; 55.68], 2, [t '2 gens, one infinite range']);

    mpc.gen(1, [QMIN QMAX]) = [-50 Inf];
    mpc.gen(2, [QMIN QMAX]) = [50 Inf];
    r = runpf(mpc, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.gen(1:2, QG), [-33.12; 57.18], 2, [t '2 gens, both infinite range']);

    mpc.gen(1, [QMIN QMAX]) = [-50 Inf];
    mpc.gen(2, [QMIN QMAX]) = [-Inf 150];
    r = runpf(mpc, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.gen(1:2, QG), [76.07; -52], 2, [t '2 gens, both infinite range']);

    mpc.gen(1, [QMIN QMAX]) = [-Inf Inf];
    mpc.gen(2, [QMIN QMAX]) = [-Inf Inf];
    r = runpf(mpc, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.gen(1:2, QG), [12.03; 12.03], 2, [t '2 gens, both infinite range']);

    t = sprintf('%s - reactive generation allocation : ', cfg{k}{1});
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
    t_is(r.gen(2:4, QG), [-5; -5; 10] + [1; 2; 3]*1.989129794, 7, [t 'PV bus']);
    t_is(r.gen(5:7, QG), [1; 2; -3], 8, [t 'PQ bus']);

    %% network with islands
    t = sprintf('%s - network w/islands : AC PF : ', cfg{k}{1});
    mpc = mpc1;
    r = runpf(mpc, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.bus( 1:9,  VA), bus_soln(:, VA), 7, [t 'voltage angles 1']);
    t_is(r.bus(10:18, VA), bus_soln(:, VA), 7, [t 'voltage angles 2']);
    Pg = [gen_soln(1, PG)-30; 30; gen_soln(2:3, PG)];
    t_is(r.gen(1:4, PG), Pg, 6, [t 'active power generation 1']);
    t_is(r.gen(5:8, PG), Pg, 6, [t 'active power generation 2']);

    t = sprintf('%s - all buses isolated : ', cfg{k}{1});
    mpc.bus(:, BUS_TYPE) = NONE;
    try
        r = runpf(mpc, mpopt);
        t_is(r.success, 0, 12, [t 'success = 0']);
    catch
        t_ok(0, [t 'unexpected fatal error']);
    end

    %% case 14 with Q limits
    t = sprintf('%s - pf.enforce_q_lims == 0 : ', cfg{k}{1});
    mpc = loadcase('case14');
    mpc.gen(1, QMIN) = -10;
    mpc.gen(:, QMAX) = [10; 30; 29; 15; 15];
    bt0 = mpc.bus(:, BUS_TYPE);
    bt = bt0;
    mpopt = mpoption(mpopt, 'pf.enforce_q_lims', 0);
    r = runpf(mpc, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.bus(:, BUS_TYPE), bt, 12, [t 'bus type']);
    t_is(r.gen(:, QG), [-16.549300542; 43.557100134; 25.075348495; 12.730944405; 17.623451366], 6, [t 'Qg']);

    t = sprintf('%s - pf.enforce_q_lims == 1 : ', cfg{k}{1});
    mpopt = mpoption(mpopt, 'pf.enforce_q_lims', 1);
    r = runpf(mpc, mpopt);
    bt = bt0;
    bt([1 2 3 8]) = [PQ PQ REF PQ];
    t_is(r.success, 0, 12, [t 'success = 0']);
    t_is(r.bus(:, BUS_TYPE), bt, 12, [t 'bus type']);
    t_is(r.gen(:, QG), [-10; 30; 31.608422873; 16.420423190; 15], 4, [t 'Qg']);

    t = sprintf('%s - pf.enforce_q_lims == 2 : ', cfg{k}{1});
    mpopt = mpoption(mpopt, 'pf.enforce_q_lims', 2);
    r = runpf(mpc, mpopt);
    bt = bt0;
    bt([1 2 3 6 8]) = [REF PQ PQ PQ PQ];
    t_ok(r.success, [t 'success']);
    t_is(r.bus(:, BUS_TYPE), bt, 12, [t 'bus type']);
    t_is(r.gen(:, QG), [-6.30936644; 30; 29; 15; 15], 5, [t 'Qg']);
  end
end

t = sprintf('%s - pf.enforce_q_lims == 1 : ', 'NR-SP');
mpopt = mpoption(mpopt, 'pf.alg', 'NR-SP', 'pf.tol', 1e-9);
mpopt = mpoption(mpopt, 'pf.enforce_q_lims', 1);
r = runpf('case14', mpopt);
t_ok(r.success, [t 'success']);
t_is(r.iterations, 5, 12, [t 'iterations']);
t_is(r.bus(1, VA), 0, 12, [t 'ref bus Va = 0']);

t = sprintf('%s - pf.enforce_q_lims == 1 : ', 'NR-IC');
mpopt = mpoption(mpopt, 'pf.alg', 'NR-SP', 'pf.tol', 1e-9);
mpopt = mpoption(mpopt, 'pf.enforce_q_lims', 1);
r = runpf('case14', mpopt);
t_ok(r.success, [t 'success']);
t_is(r.iterations, 5, 12, [t 'iterations']);
t_is(r.bus(1, VA), 0, 12, [t 'ref bus Va = 0']);

t_end;

if have_feature('octave')
    warning(s1.state, file_in_path_warn_id);
end
