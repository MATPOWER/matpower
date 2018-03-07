function t_opf_softlims2(quiet)
%T_OPF_SOFTLIMS  Tests for userfcn callbacks (softlims) w/OPF.
%   Includes high-level tests of soft limits implementations.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   updated by Eran Schweitzer, Eran.Schweitzer@asu.edu
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

casefile = 'case9';
if quiet
    verbose = 0;
else
    verbose = 0;
end

% t_begin(59+37*4, quiet);
t_begin(629, quiet);
%% define constants
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

%% set options
mpopt = mpoption('opf.dc.solver', 'MIPS', 'opf.ac.solver', 'IPOPT');
% mpopt = mpoption('opf.dc.solver', 'MIPS', 'opf.ac.solver', 'MIPS');
%mpopt = mpoption('opf.dc.solver', 'GLPK', 'opf.ac.solver', 'FMINCON');
%mpopt = mpoption('opf.dc.solver', 'GUROBI', 'opf.ac.solver', 'IPOPT');
%mpopt = mpoption('opf.dc.solver', 'CPLEX', 'opf.ac.solver', 'KNITRO');
%mpopt = mpoption(mpopt, 'fmincon.tol_x', 1e-9, 'fmincon.tol_f', 1e-9);
% mpopt = mpoption(mpopt, 'ipopt.opts.tol', 1e-10, 'ipopt.opts.acceptable_tol', 1e-10);
%mpopt = mpoption(mpopt, 'knitro.tol_x', 1e-9, 'knitro.tol_f', 1e-9);
mpopt = mpoption(mpopt, 'opf.violation', 1e-6, 'mips.gradtol', 1e-8, ...
        'mips.comptol', 1e-8, 'mips.costtol', 1e-9);
if verbose <= 1
    mpopt = mpoption(mpopt, 'out.all', 0);
end
mpopt = mpoption(mpopt, 'verbose', verbose);

if have_fcn('matlab')
    sing_matrix_warn_id = 'MATLAB:nearlySingularMatrix';
    s2 = warning('query', sing_matrix_warn_id);
    warning('off', sing_matrix_warn_id);
end

%% load and modify case file
mpc = loadcase(casefile);
nl = size(mpc.branch, 1);
mpc.gencost(:, NCOST) = 2;
mpc.gencost(:, COST) = [];
mpc.gencost(:, COST) = [50; 40; 25];
mpc.gencost(:, COST+1) = 0;
mpc.branch(2:6, RATE_A) = [0; 120; 0; 120; 0];
mpc.branch(3, RATE_A) = 120;
mpc.branch(5, RATE_A) = 120;
mpc.branch = [mpc.branch(1:4, :); mpc.branch(4, :);  mpc.branch(5:end, :)];
mpc.branch(5, BR_STATUS) = 0;

mpc.branch(:, ANGMIN) = -60;
mpc.branch(:, ANGMAX) = 60;

nl = size(mpc.branch, 1);   %% number of branches

%% create soft limit inputs (emulates previous implementation just for branches)
mpc.softlims.RATE_A.type = 'unbnd';
mpc.softlims.RATE_A.idx  = (2:nl)';
mpc.softlims.RATE_A.cost = 100 * ones(nl-1, 1);

for prop = {'VMIN', 'VMAX', 'ANGMIN', 'ANGMAX', 'PMIN', 'PMAX', 'QMIN', 'QMAX'}
    prop = prop{1};
    mpc.softlims.(prop).type = 'none';
end
mpc0 = mpc;     %% save to initialize later runs

%% Default softlims structrue
sdefault = struct();
for prop = {'RATE_A','VMIN', 'VMAX', 'ANGMIN', 'ANGMAX', 'PMIN', 'PMAX', 'QMIN', 'QMAX'}
    prop = prop{1};
    if ismember(prop, {'VMAX', 'VMIN'})
        sdefault.(prop).type = 'frac';
        sdefault.(prop).ub     = 0.5;
    elseif ismember(prop, {'ANGMAX', 'ANGMIN'})
        sdefault.(prop).type = 'cnst';
        sdefault.(prop).ub   = 360;
    elseif strcmp(prop, 'RATE_A')
        sdefault.(prop).type = 'frac';
        sdefault.(prop).ub   = 0.5;
    elseif ismember(prop, {'PMAX', 'QMAX', 'QMIN'})
        sdefault.(prop).type = 'unbnd';
    elseif strcmp(prop, 'PMIN')
        sdefault.(prop).type = 'frac';
        sdefault.(prop).ub   = 1;
    end
end
%% Default settings check
mpc = mpc0;
mpc.softlims = sdefault;

t = 'Defaults test - no active softlims: ';
r = toggle_run_check(mpc, mpopt, t, 0);
delta = calc_branch_angle(r);
overload_loop(r, t, 0);
t_is(r.bus(:,VM), [1.040963;1.100000;1.099999;1.042513;1.041787;1.096141;1.079794;1.089516;1.031234], 4, [t 'Vm']);
t_is(r.branch(:,PF), [17.3121  -25.0634 -115.1668  239.9012         0  119.9646   18.5301  -66.4562   84.9411  -42.2020].', 4, [t 'PF'])
t_is(r.branch(:,QF), [-2.7212   -2.8421  -16.2418   21.2760         0   -2.9141  -25.3254  -17.1234    8.9449  -17.4049].', 4, [t 'QF'])
t_is(r.branch(:,PT), [-17.3121   25.1668  119.9366 -239.9012         0 -118.5301  -18.4848   66.4562  -82.7980   42.3756].',4, [t 'PT'])
t_is(r.branch(:,QT), [2.8844  -13.7582   -3.9017    6.8158         0   -9.6746    8.1785   19.6031  -32.5951   -0.0423].', 4, [t 'QT'])
t_is(delta, [0.5265   -1.2681   -9.9352    6.6955         0    5.8081    0.7187   -1.9861    6.5458   -1.8692].', 4, [t 'delta'])
t_is(r.gen(:,PG), [17.3121   66.4562  239.9012].', 4, [t 'PG'])
t_is(r.gen(:,QG), [-2.7212   19.6031   21.2760].', 4, [t 'QG'])
t_is(r.bus(:,MU_VMAX), [0.000000;526.455671;232.322240;0.000000;0.000000;0.000000;0.000000;0.000000;0.000000], 4, [t 'mu VM ub'])
t_is(r.bus(:,MU_VMIN), zeros(9,1), 4, [t 'mu VM lb'])
t_is(sum(r.branch(:,MU_SF:MU_ST),2), [0.0000         0   29.2428         0         0    8.5610         0    0.0000    0.0000    0.0000].', 4, [t 'mu SF+ST'])
t_is(r.branch(:,MU_ANGMAX), zeros(10,1), 4, [t 'mu ANGMAX'])
t_is(r.branch(:,MU_ANGMIN), zeros(10,1), 4, [t 'mu ANGMIN'])
t_is(r.gen(:,MU_PMAX), zeros(3,1), 4, [t 'mu PMAX'])
t_is(r.gen(:,MU_PMIN), zeros(3,1), 4, [t 'mu PMIN'])
t_is(r.gen(:,MU_QMAX), zeros(3,1), 4, [t 'mu QMAX'])
t_is(r.gen(:,MU_QMIN), zeros(3,1), 4, [t 'mu QMIN'])

t = 'Defaults test - softlims on: ';
mpc = toggle_softlims(mpc,'on');
r = toggle_run_check(mpc, mpopt, t, 1);
delta = calc_branch_angle(r);
overload_loop(r, t, 1);
t_is(r.bus(:,VM), [1.040963;1.100000;1.099999;1.042513;1.041787;1.096141;1.079794;1.089516;1.031234], 4, [t 'Vm']);
t_is(r.branch(:,PF), [17.3121  -25.0634 -115.1668  239.9012         0  119.9646   18.5301  -66.4562   84.9411  -42.2020].', 4, [t 'PF'])
t_is(r.branch(:,QF), [-2.7212   -2.8421  -16.2418   21.2760         0   -2.9141  -25.3254  -17.1234    8.9449  -17.4049].', 4, [t 'QF'])
t_is(r.branch(:,PT), [-17.3121   25.1668  119.9366 -239.9012         0 -118.5301  -18.4848   66.4562  -82.7980   42.3756].',4, [t 'PT'])
t_is(r.branch(:,QT), [2.8844  -13.7582   -3.9017    6.8158         0   -9.6746    8.1785   19.6031  -32.5951   -0.0423].', 4, [t 'QT'])
t_is(delta, [0.5265   -1.2681   -9.9352    6.6955         0    5.8081    0.7187   -1.9861    6.5458   -1.8692].', 4, [t 'delta'])
t_is(r.gen(:,PG), [17.3121   66.4562  239.9012].', 4, [t 'PG'])
t_is(r.gen(:,QG), [-2.7212   19.6031   21.2760].', 4, [t 'QG'])
t_is(r.bus(:,MU_VMAX), [0.000000;526.455671;232.322240;0.000000;0.000000;0.000000;0.000000;0.000000;0.000000], 4, [t 'mu VM ub'])
t_is(r.bus(:,MU_VMIN), zeros(9,1), 4, [t 'mu VM lb'])
t_is(sum(r.branch(:,MU_SF:MU_ST),2), [0.0000         0   29.2428         0         0    8.5610         0    0.0000    0.0000    0.0000].', 4, [t 'mu SF+ST'])
t_is(r.branch(:,MU_ANGMAX), zeros(10,1), 4, [t 'mu ANGMAX'])
t_is(r.branch(:,MU_ANGMIN), zeros(10,1), 4, [t 'mu ANGMIN'])
t_is(r.gen(:,MU_PMAX), zeros(3,1), 4, [t 'mu PMAX'])
t_is(r.gen(:,MU_PMIN), zeros(3,1), 4, [t 'mu PMIN'])
t_is(r.gen(:,MU_QMAX), zeros(3,1), 4, [t 'mu QMAX'])
t_is(r.gen(:,MU_QMIN), zeros(3,1), 4, [t 'mu QMIN'])

%%% tighten limits to get violations
t = 'Defaults test - softlimits with overloads: ';
mpc = mpc0;
mpc.softlims = sdefault;
mpc.bus(4,VMAX) = 0.9;
mpc.softlims.VMAX.cost = 75000;
mpc.bus(9,VMIN)  = 1.05;
mpc.softlims.VMIN.cost = 75000;
mpc.softlims.RATE_A.cost = 10;
mpc.gen(1,PMIN) = 25;
mpc.softlims.PMIN.cost = 5;
mpc.gen(3,PMAX) = 225;
mpc.softlims.PMAX.idx  = [1;2;3];
mpc.softlims.PMAX.cost = [1000;1000;3];
mpc.gen(2,QMAX) = 50;
mpc.softlims.QMAX.idx  = [1;2;3];
mpc.softlims.QMAX.cost = [1000; 20; 1000];
mpc.gen(1,QMIN) = -45;
mpc.softlims.QMIN.idx  = [1;2;3];
mpc.softlims.QMIN.cost = [10; 1000; 1000];
mpc.branch(3,ANGMIN) = -5;
mpc.branch(4,ANGMAX) = +5;
mpc.softlims.ANGMIN.idx = (1:nl).';
mpc.softlims.ANGMAX.idx = (1:nl).';
mpc.softlims.ANGMIN.cost = 1e3*ones(10,1);
mpc.softlims.ANGMAX.cost = 1e3*ones(10,1);
mpc.softlims.ANGMIN.cost(3) = 20;
mpc.softlims.ANGMAX.cost(4) = 20;
mpc = toggle_softlims(mpc,'on');
r = toggle_run_check(mpc, mpopt, t, 1);
overload_loop(r, t, 1);
mu_cost_test(r,t);

% test save case
t = 'savecase(fname, r) : ';
fn = sprintf('softlims_savecase_test_%d', fix(1e8*rand));
savecase(fn, r);
mpc1 = loadcase(fn);
delete([fn '.m']);
t_ok(isfield(mpc1, 'softlims'), [t 'mpc.softlims']);
for prop = fieldnames(r.softlims).'
    t_ok(isfield(mpc1.softlims, prop{1}), [t 'mpc.softlims.' prop{1}])
    for field = fieldnames(r.softlims.(prop{1})).'
        if any(isinf(r.softlims.(prop{1}).(field{1})))
            t_ok(isequal(mpc1.softlims.(prop{1}).(field{1}), r.softlims.(prop{1}).(field{1})), [t 'mpc.softlims.' prop{1} '.' field{1}])
        else
            t_is(mpc1.softlims.(prop{1}).(field{1}), r.softlims.(prop{1}).(field{1}), 5, [t 'mpc.softlims.' prop{1} '.' field{1}])
        end
    end
end

%% defaults with DC
mpc = mpc0;
mpc.softlims = sdefault;
t = 'DC - hard limits : ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
r = rundcopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
delta = calc_branch_angle(r);
overload_loop(r, t, 0);
t_is(r.gen(:, PG), [12.687; 62.3130; 240], 4, [t 'Pg']);
t_is(r.branch(:, PF), [12.687; -30; -120; 240; 0; 120; 20; -62.3130; 82.3130; -42.687], 4, [t 'Pf']);
t_is(delta, [0.4187   -1.5814  -11.6883    8.0581         0    6.9305    0.8251   -2.2314    7.5931   -2.0789].', 4, [t 'delta'])
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 35.6504; 0; 0; 7.9756; 0; 0; 0; 0], 4, [t 'mu Pf']);


t = 'DC - softlims on: ';
mpc = toggle_softlims(mpc,'on');
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = rundcopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
delta = calc_branch_angle(r);
overload_loop(r, t, 1);
t_is(r.gen(:, PG), [12.687; 62.3130; 240], 4, [t 'Pg']);
t_is(r.branch(:, PF), [12.687; -30; -120; 240; 0; 120; 20; -62.3130; 82.3130; -42.687], 4, [t 'Pf']);
t_is(delta, [0.4187   -1.5814  -11.6883    8.0581         0    6.9305    0.8251   -2.2314    7.5931   -2.0789].', 4, [t 'delta'])
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 35.6504; 0; 0; 7.9756; 0; 0; 0; 0], 4, [t 'mu Pf']);

%%% tighten limits to get violations
t = 'DC - softlimits with overloads: ';
mpc = mpc0;
mpc.softlims = sdefault;
mpc.branch([3,6,7,8],RATE_A) = 75;
mpc.softlims.RATE_A.cost = 15;
mpc.gen(1,PMIN) = 75;
mpc.softlims.PMIN.idx = [1;2;3];
mpc.softlims.PMIN.cost = [5; 1000; 1000];
mpc.gen(3,PMAX) = 90;
mpc.softlims.PMAX.idx = [1;2;3];
mpc.softlims.PMAX.cost = [1000; 1000; 10];
mpc.branch(3,ANGMIN) = -5;
mpc.softlims.ANGMIN.idx = (1:10)';
mpc.softlims.ANGMIN.cost = 1000*ones(10,1);
mpc.softlims.ANGMIN.cost(3) = 5;
mpc.branch(9,ANGMAX) = +5;
mpc.softlims.ANGMAX.idx = (1:10)';
mpc.softlims.ANGMAX.cost = 1000*ones(10,1);
mpc.softlims.ANGMAX.cost(9) = 5;
mpc = toggle_softlims(mpc,'on');
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = rundcopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
overload_loop(r, t, 1);
mu_cost_test(r,t);
%% voltage magnitude slack
% unbounded limits
mpc = mpc0;
mpc.softlims.RATE_A.type = 'none';
mpc.softlims.VMAX.type   = 'unbnd';
mpc.softlims.VMIN.type   = 'unbnd';
t = 'Voltage limits - hard limits: ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
overload_loop(r, t, 0);
t_is(r.bus(:,VM), [1.040963;1.100000;1.099999;1.042513;1.041787;1.096141;1.079794;1.089516;1.031234], 4, [t 'Vm']);
t_is(r.bus(:,MU_VMAX), [0.000000;526.455671;232.322240;0.000000;0.000000;0.000000;0.000000;0.000000;0.000000], 4, [t 'mu VM ub'])

mpc = toggle_softlims(mpc,'on');
t = 'Voltage limits (unbounded) - soft limits (satisfied): ';
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.bus(:,VM), [1.040963;1.100000;1.099999;1.042513;1.041787;1.096141;1.079794;1.089516;1.031234], 4, [t 'Vm']);
t_is(r.bus(:,MU_VMAX), [0.000000;526.455671;232.322240;0.000000;0.000000;0.000000;0.000000;0.000000;0.000000], 4, [t 'mu VM ub'])


mpc.bus(4:9, [VMIN, VMAX]) = 1;
t = 'voltage limits (unbounded) - softlimits with overloads: ';
r = toggle_run_check(mpc, mpopt, t, 1);
overload_loop(r, t, 1, {'VMIN', 'VMAX'})
overload_loop(r, t, 0, {'VMIN', 'VMAX'})
mu_cost_test(r,t);


% fractional limits
mpc = mpc0;
mpc.softlims = sdefault;
for prop = fieldnames(mpc.softlims).'
    if ~ismember(prop{1},{'VMIN','VMAX'})
        mpc.softlims.(prop{1}).type = 'none';
    end
end
mpc.bus(4:9, [VMIN, VMAX]) = 1;
mpc = toggle_softlims(mpc,'on');
t = 'voltage limits (fractional) - softlimits with overloads: ';
r = toggle_run_check(mpc, mpopt, t, 1);
overload_loop(r, t, 1, {'VMIN', 'VMAX'})
overload_loop(r, t, 0, {'VMIN', 'VMAX'})
mu_cost_test(r,t);

% constant upper bound limit
mpc = mpc0;
mpc.softlims = sdefault;
for prop = fieldnames(mpc.softlims).'
    if strcmp(prop{1},'VMIN')
        mpc.softlims.VMIN.type = 'cnst';
        mpc.softlims.VMIN.ub   = 0.7;
    elseif strcmp(prop{1},'VMAX')
        mpc.softlims.VMAX.type = 'cnst';
        mpc.softlims.VMAX.ub   = 1.2;
    else
        mpc.softlims.(prop{1}).type = 'none';
    end
end
mpc.bus(4:9, [VMIN, VMAX]) = 1;
mpc = toggle_softlims(mpc,'on');
t = 'voltage limits (constant) - softlimits with overloads: ';
r = toggle_run_check(mpc, mpopt, t, 1);
overload_loop(r, t, 1, {'VMIN', 'VMAX'})
overload_loop(r, t, 0, {'VMIN', 'VMAX'})
mu_cost_test(r,t);

%% Tests
% the following all all the test from the original softlims implementation
% that only included RATE_A. The only changes are adaptations to the
% changed data structure.
mpc = mpc0;
t = 'DC - hard limits : ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
r = rundcopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(~isfield(r.softlims.RATE_A, 'overload'), [t 'no softlims.RATE_A.overload']);
t_ok(~isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'no softlims.RATE_A.ovl_cost']);
t_is(r.f, 9126.87, 4, [t 'f']);
t_is(r.gen(:, PG), [12.687; 62.3130; 240], 4, [t 'Pg']);
t_is(r.branch(:, PF), [12.687; -30; -120; 240; 0; 120; 20; -62.3130; 82.3130; -42.687], 4, [t 'Pf']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 35.6504; 0; 0; 7.9756; 0; 0; 0; 0], 4, [t 'mu Pf']);

t = 'DC - soft limits (satisfied) : ';
mpc = toggle_softlims(mpc, 'on');
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = rundcopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9126.87, 4, [t 'f']);
t_is(r.gen(:, PG), [12.687; 62.3130; 240], 4, [t 'Pg']);
t_is(r.branch(:, PF), [12.687; -30; -120; 240; 0; 120; 20; -62.3130; 82.3130; -42.687], 4, [t 'Pf']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 35.6504; 0; 0; 7.9756; 0; 0; 0; 0], 4, [t 'mu Pf']);
t_is(r.softlims.RATE_A.overload, zeros(10, 1), 12, [t 'softlims.overload']);
t_is(r.softlims.RATE_A.ovl_cost, zeros(10, 1), 12, [t 'softlims.ovl_cost']);
t_is(r.order.branch.status.on(r.order.int.softlims.RATE_A.idx), [3; 6; 8; 9; 10], 12, [t 'mu Pf']);

t = 'savecase(fname, mpc) : ';
fn = sprintf('softlims_savecase_test_%d', fix(1e8*rand));
savecase(fn, mpc);
mpc1 = loadcase(fn);
delete([fn '.m']);
t_ok(isfield(mpc1, 'softlims'), [t 'mpc.softlims']);
for prop = fieldnames(mpc.softlims).'
    t_ok(isfield(mpc1.softlims, prop{1}), [t 'mpc.softlims.' prop{1}])
    for field = fieldnames(mpc.softlims.(prop{1})).'
        if any(isinf(mpc.softlims.(prop{1}).(field{1})))
            t_ok(isequal(mpc1.softlims.(prop{1}).(field{1}), mpc.softlims.(prop{1}).(field{1})), [t 'mpc.softlims.' prop{1} '.' field{1}])
        else
            t_is(mpc1.softlims.(prop{1}).(field{1}), mpc.softlims.(prop{1}).(field{1}), 5, [t 'mpc.softlims.' prop{1} '.' field{1}])
        end
    end
end
 
t = 'savecase(fname, results) + gentype/genfuel: ';
r.genfuel = {'ng'; 'coal'; 'hydro'};
r.gentype = {'GT'; 'ST'; 'HY'};
fn = sprintf('softlims_savecase_test_%d', fix(1e8*rand));
savecase(fn, r);
mpc1 = loadcase(fn);
delete([fn '.m']);
t_ok(isfield(mpc1, 'softlims'), [t 'results.softlims']);
for prop = fieldnames(r.softlims).'
    t_ok(isfield(mpc1.softlims, prop{1}), [t 'results.softlims.' prop{1}])
    for field = fieldnames(r.softlims.(prop{1})).'
        if any(isinf(r.softlims.(prop{1}).(field{1})))
            t_ok(isequal(mpc1.softlims.(prop{1}).(field{1}), r.softlims.(prop{1}).(field{1})), [t 'results.softlims.' prop{1} '.' field{1}])
        else
            t_is(mpc1.softlims.(prop{1}).(field{1}), r.softlims.(prop{1}).(field{1}), 5, [t 'results.softlims.' prop{1} '.' field{1}])    
        end
    end
end
t_ok(isfield(mpc1, 'gentype'), [t 'results.gentype']);
t_ok(isfield(mpc1, 'genfuel'), [t 'results.genfuel']);
if isfield(mpc1, 'gentype')
    t_ok(isequal(mpc1.gentype, r.gentype), [t 'results.gentype']);
else
    t_ok(0, [t 'results.gentype']);
end
if isfield(mpc1, 'genfuel')
    t_ok(isequal(mpc1.genfuel, r.genfuel), [t 'results.genfuel']);
else
    t_ok(0, [t 'results.genfuel']);
end

t = 'DC - soft limits (violated) : ';
mpc = rmfield(mpc, 'softlims');
mpc.softlims.RATE_A.type = 'unbnd';
mpc.softlims.RATE_A.cost = 20;
for prop = {'VMIN', 'VMAX', 'ANGMIN', 'ANGMAX', 'PMIN', 'PMAX', 'QMIN', 'QMAX'}
    prop = prop{1};
    mpc.softlims.(prop).type = 'none';
end
% mpc.softlims.cost = 20;
r = rundcopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9106.5059, 4, [t 'f']);
t_is(r.gen(:, PG), [10; 63.6988; 241.3012], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -31.3012; -121.3012; 241.3012; 0; 120; 20; -63.6988; 83.6988; -41.3012], 4, [t 'Pf']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 20; 0; 0; 13.2992; 0; 0; 0; 0], 4, [t 'mu Pf']);
t_is(r.softlims.RATE_A.overload, [0; 0; 1.3011811; 0; 0; 0; 0; 0; 0; 0], 6, [t 'softlims.RATE_A.overload']);
t_is(r.softlims.RATE_A.ovl_cost, [0; 0; 26.023622; 0; 0; 0; 0; 0; 0; 0], 6, [t 'softlims.RATE_A.ovl_cost']);

t = 'savecase(fname, mpc) : ';
fn = sprintf('softlims_savecase_test_%d', fix(1e8*rand));
savecase(fn, mpc);
mpc1 = loadcase(fn);
delete([fn '.m']);
t_ok(isfield(mpc1, 'softlims'), [t 'mpc.softlims']);
for prop = fieldnames(mpc.softlims).'
    t_ok(isfield(mpc1.softlims, prop{1}), [t 'mpc.softlims.' prop{1}])
    for field = fieldnames(mpc.softlims.(prop{1})).'
        if any(isinf(mpc.softlims.(prop{1}).(field{1})))
            t_ok(isequal(mpc1.softlims.(prop{1}).(field{1}), mpc.softlims.(prop{1}).(field{1})), [t 'mpc.softlims.' prop{1} '.' field{1}])
        else
            t_is(mpc1.softlims.(prop{1}).(field{1}), mpc.softlims.(prop{1}).(field{1}), 5, [t 'mpc.softlims.' prop{1} '.' field{1}])
        end
    end
end

t = 'savecase(fname, results) : ';
fn = sprintf('softlims_savecase_test_%d', fix(1e8*rand));
savecase(fn, r);
mpc1 = loadcase(fn);
delete([fn '.m']);
t_ok(isfield(mpc1, 'softlims'), [t 'results.softlims']);
for prop = fieldnames(r.softlims).'
    t_ok(isfield(mpc1.softlims, prop{1}), [t 'results.softlims.' prop{1}])
    for field = fieldnames(r.softlims.(prop{1})).'
        if any(isinf(r.softlims.(prop{1}).(field{1})))
            t_ok(isequal(mpc1.softlims.(prop{1}).(field{1}), r.softlims.(prop{1}).(field{1})), [t 'results.softlims.' prop{1} '.' field{1}])
        else
            t_is(mpc1.softlims.(prop{1}).(field{1}), r.softlims.(prop{1}).(field{1}), 5, [t 'results.softlims.' prop{1} '.' field{1}])   
        end
    end
end

%% -----  AC OPF (opf.flow_lim = 'S' -----
mpc = mpc0;
t = 'AC (flow_lim ''S'') - hard limits : ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(~isfield(r.softlims.RATE_A, 'overload'), [t 'no softlims.RATE_A.overload']);
t_ok(~isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'no softlims.RATE_A.ovl_cost']);
t_is(r.f, 9521.385584, 4, [t 'f']);
t_is(r.gen(:, PG), [17.312141; 66.456235; 239.901164], 4, [t 'Pg']);
t_is(r.branch(:, PF), [17.312141; -25.063413; -115.166831; 239.901164; 0; 119.964612; 18.530059; -66.456235; 84.941080; -42.201990], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-17.312141; 25.166831; 119.936552; -239.901164; 0; -118.530059; -18.484844; 66.456235; -82.798010; 42.375554], 4, [t 'Pt']);
t_is(r.branch(:, QF), [-2.721149; -2.842148; -16.241771; 21.275974; 0; -2.914100; -25.325409; -17.123405; 8.944864; -17.404939], 4, [t 'Qf']);
t_is(r.branch(:, QT), [2.884399; -13.758229; -3.901717; 6.815818; 0; -9.674591; 8.178541; 19.603112; -32.595061; -0.042250], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 29.242772; 0; 0; 8.561015; 0; 0; 0; 0], 4, [t 'mu Sf+St']);

t = 'AC (flow_lim ''S'') - soft limits (satisfied) : ';
mpc = toggle_softlims(mpc, 'on');
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9521.385584, 4, [t 'f']);
t_is(r.gen(:, PG), [17.312141; 66.456235; 239.901164], 4, [t 'Pg']);
t_is(r.branch(:, PF), [17.312141; -25.063413; -115.166831; 239.901164; 0; 119.964612; 18.530059; -66.456235; 84.941080; -42.201990], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-17.312141; 25.166831; 119.936552; -239.901164; 0; -118.530059; -18.484844; 66.456235; -82.798010; 42.375554], 4, [t 'Pt']);
t_is(r.branch(:, QF), [-2.721149; -2.842148; -16.241771; 21.275974; 0; -2.914100; -25.325409; -17.123405; 8.944864; -17.404939], 4, [t 'Qf']);
t_is(r.branch(:, QT), [2.884399; -13.758229; -3.901717; 6.815818; 0; -9.674591; 8.178541; 19.603112; -32.595061; -0.042250], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 29.242772; 0; 0; 8.561015; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.RATE_A.overload, zeros(10, 1), 12, [t 'softlims.overload']);
t_is(r.softlims.RATE_A.ovl_cost, zeros(10, 1), 12, [t 'softlims.ovl_cost']);
t_is(r.order.branch.status.on(r.order.int.softlims.RATE_A.idx), [3; 6; 8; 9; 10], 12, [t 'softlims.RATE_A.idx']);

t = 'AC (flow_lim ''S'') - soft limits (violated) : ';
mpc = rmfield(mpc, 'softlims');
mpc.softlims.RATE_A.type = 'unbnd';
mpc.softlims.RATE_A.cost = 20;
for prop = {'VMIN', 'VMAX', 'ANGMIN', 'ANGMAX', 'PMIN', 'PMAX', 'QMIN', 'QMAX'}
    prop = prop{1};
    mpc.softlims.(prop).type = 'none';
end
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9486.265634, 4, [t 'f']);
t_is(r.gen(:, PG), [10; 70.325987; 243.758931], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -28.621625; -118.761245; 243.758931; 0; 119.958848; 18.527937; -70.325987; 88.808549; -38.472172], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-10; 28.761245; 123.8; -243.758931; 0; -118.527937; -18.482561; 70.325987; -86.527828; 38.621625], 4, [t 'Pt']);
t_is(r.branch(:, QF), [3.364131; 0.552145; -12.847555; 19.471020; 0; -3.142423; -25.464124; -14.475416; 6.202740; -20.630111], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-3.306136; -17.152445; -6.346356; 9.488779; 0; -9.535876; 8.272676; 17.182536; -29.369889; 2.753991], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 20; 0; 0; 11.575125; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.RATE_A.overload, [0; 0; 3.962643; 0; 0; 0; 0; 0; 0; 0], 6, [t 'softlims.RATE_A.overload']);
t_is(r.softlims.RATE_A.ovl_cost, [0; 0; 79.252860; 0; 0; 0; 0; 0; 0; 0], 4, [t 'softlims.RATE_A.ovl_cost']);


%% -----  AC OPF (opf.flow_lim = 'I' -----
mpopt = mpoption(mpopt, 'opf.flow_lim', 'I');
mpc = mpc0;
t = 'AC (flow_lim ''I'') - hard limits : ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(~isfield(r.softlims.RATE_A, 'overload'), [t 'no softlims.RATE_A.overload']);
t_ok(~isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'no softlims.RATE_A.ovl_cost']);
t_is(r.f, 9178.300537, 4, [t 'f']);
t_is(r.gen(:, PG), [10; 54.432302; 260.040338], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -32.676659; -122.910718; 260.040338; 0; 131.831065; 30.115930; -54.432302; 84.451894; -42.475047], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-10; 32.910718; 128.209273; -260.040338; 0; -130.115930; -30.019591; 54.432302; -82.524953; 42.676659], 4, [t 'Pt']);
t_is(r.branch(:, QF), [27.824349; 14.237742; 1.341441; 7.733221; 0; -4.887428; -29.462755; -4.743706; -7.797128; -31.778587], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-27.408203; -31.341441; -20.438656; 25.326084; 0; -5.537245; 12.540834; 6.294583; -18.221413; 13.170461], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 0; 0; 0; 20.114712; 0; 0; 0; 0], 4, [t 'mu Sf+St']);

t = 'AC (flow_lim ''I'') - soft limits (satisfied) : ';
mpc = toggle_softlims(mpc, 'on');
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9178.300537, 4, [t 'f']);
t_is(r.gen(:, PG), [10; 54.432302; 260.040338], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -32.676659; -122.910718; 260.040338; 0; 131.831065; 30.115930; -54.432302; 84.451894; -42.475047], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-10; 32.910718; 128.209273; -260.040338; 0; -130.115930; -30.019591; 54.432302; -82.524953; 42.676659], 4, [t 'Pt']);
t_is(r.branch(:, QF), [27.824349; 14.237742; 1.341441; 7.733221; 0; -4.887428; -29.462755; -4.743706; -7.797128; -31.778587], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-27.408203; -31.341441; -20.438656; 25.326084; 0; -5.537245; 12.540834; 6.294583; -18.221413; 13.170461], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 0; 0; 0; 20.114712; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.RATE_A.overload, zeros(10, 1), 12, [t 'softlims.RATE_A.overload']);
t_is(r.softlims.RATE_A.ovl_cost, zeros(10, 1), 12, [t 'softlims.RATE_A.ovl_cost']);
t_is(r.order.branch.status.on(r.order.int.softlims.RATE_A.idx), [3; 6; 8; 9; 10], 12, [t 'softlims.RATE_A.idx']);

t = 'AC (flow_lim ''I'') - soft limits (violated) : ';
mpc = rmfield(mpc, 'softlims');
mpc.softlims.RATE_A.type = 'unbnd';
mpc.softlims.RATE_A.cost = 15;
for prop = {'VMIN', 'VMAX', 'ANGMIN', 'ANGMAX', 'PMIN', 'PMAX', 'QMIN', 'QMAX'}
    prop = prop{1};
    mpc.softlims.(prop).type = 'none';
end
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9146.413838, 4, [t 'f']);
t_is(r.gen(:, PG), [10; 44.862996; 270], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -34.955827; -125.208395; 270; 0; 139.280487; 37.364683; -44.862996; 82.094229; -44.742338], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-10; 35.208395; 130.719513; -270; 0; -137.364683; -37.231233; 44.862996; -80.257662; 44.955827], 4, [t 'Pt']);
t_is(r.branch(:, QF), [25.262645; 13.437313; 0.303801; 13.552057; 0; -3.645909; -29.947558; -7.234346; -6.145528; -29.826268], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-24.907470; -30.303801; -18.341454; 21.987363; 0; -5.052442; 13.379873; 8.309624; -20.173732; 11.470157], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 6.260861; 0; 0; 15; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.RATE_A.overload, [0; 0; 0; 0; 0; 6.792934; 0; 0; 0; 0], 6, [t 'softlims.RATE_A.overload']);
t_is(r.softlims.RATE_A.ovl_cost, [0; 0; 0; 0; 0; 101.894009; 0; 0; 0; 0], 4, [t 'softlims.RATE_A.ovl_cost']);

%%-----  AC OPF (opf.flow_lim = '2' -----
mpopt = mpoption(mpopt, 'opf.flow_lim', '2');
mpc = mpc0;
t = 'AC (flow_lim ''2'') - hard limits : ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(~isfield(r.softlims.RATE_A, 'overload'), [t 'no softlims.RATE_A.overload']);
t_ok(~isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'no softlims.RATE_A.ovl_cost']);
t_is(r.f, 9513.613051, 4, [t 'f']);
t_is(r.gen(:, PG), [17.759246; 65.641269; 240], 4, [t 'Pg']);
t_is(r.branch(:, PF), [17.759246; -25.216499; -115.346644; 240; 0; 120; 18.577311; -65.641269; 84.170657; -42.784248], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-17.759246; 25.346644; 120; -240; 0; -118.577311; -18.529387; 65.641269; -82.215752; 42.975745], 4, [t 'Pt']);
t_is(r.branch(:, QF), [15.646152; 6.577259; -6.087866; 8.021258; 0; -4.520662; -26.623563; -6.687069; -2.629290; -27.085546], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-15.370220; -23.912134; -15.548682; 20.069344; 0; -8.376437; 9.316359; 8.954090; -22.914454; 8.792961], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 29.248033; 0; 0; 8.417115; 0; 0; 0; 0], 4, [t 'mu Sf+St']);

t = 'AC (flow_lim ''2'') - soft limits (satisfied) : ';
mpc = toggle_softlims(mpc, 'on');
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9513.613051, 4, [t 'f']);
t_is(r.gen(:, PG), [17.759246; 65.641269; 240], 4, [t 'Pg']);
t_is(r.branch(:, PF), [17.759246; -25.216499; -115.346644; 240; 0; 120; 18.577311; -65.641269; 84.170657; -42.784248], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-17.759246; 25.346644; 120; -240; 0; -118.577311; -18.529387; 65.641269; -82.215752; 42.975745], 4, [t 'Pt']);
t_is(r.branch(:, QF), [15.646152; 6.577259; -6.087866; 8.021258; 0; -4.520662; -26.623563; -6.687069; -2.629290; -27.085546], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-15.370220; -23.912134; -15.548682; 20.069344; 0; -8.376437; 9.316359; 8.954090; -22.914454; 8.792961], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 29.248033; 0; 0; 8.417115; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.RATE_A.overload, zeros(10, 1), 12, [t 'softlims.RATE_A.overload']);
t_is(r.softlims.RATE_A.ovl_cost, zeros(10, 1), 12, [t 'softlims.RATE_A.ovl_cost']);
t_is(r.order.branch.status.on(r.order.int.softlims.RATE_A.idx), [3; 6; 8; 9; 10], 12, [t 'softlims.RATE_A.idx']);

t = 'AC (flow_lim ''2'') - soft limits (violated) : ';
mpc = rmfield(mpc, 'softlims');
mpc.softlims.RATE_A.type = 'unbnd';
mpc.softlims.RATE_A.cost = 20;
for prop = {'VMIN', 'VMAX', 'ANGMIN', 'ANGMAX', 'PMIN', 'PMAX', 'QMIN', 'QMAX'}
    prop = prop{1};
    mpc.softlims.(prop).type = 'none';
end
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9476.854715, 4, [t 'f']);
t_is(r.gen(:, PG), [10; 69.833910; 244.077740; ], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -28.933383; -119.111239; 244.077740; 0; 120; 18.578292; -69.833910; 88.362736; -38.762267], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-10; 29.111239; 124.077740; -244.077740; 0; -118.578292; -18.528826; 69.833910; -86.237733; 38.933383], 4, [t 'Pt']);
t_is(r.branch(:, QF), [22.204979; 10.325903; -2.434897; 5.980065; 0; -5.159694; -27.240912; -4.815545; -5.105567; -30.250169], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-21.917920; -27.565103; -17.970864; 23.130558; 0; -7.759088; 9.921112; 7.362539; -19.749831; 11.592017], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 20; 0; 0; 11.526783; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.RATE_A.overload, [0; 0; 4.07774; 0; 0; 0; 0; 0; 0; 0], 6, [t 'softlims.RATE_A.overload']);
t_is(r.softlims.RATE_A.ovl_cost, [0; 0; 81.554809; 0; 0; 0; 0; 0; 0; 0], 4, [t 'softlims.RATE_A.ovl_cost']);

%%-----  AC OPF (opf.flow_lim = 'P' -----
mpopt = mpoption(mpopt, 'opf.flow_lim', 'P');
mpc = mpc0;
t = 'AC (flow_lim ''P'') - hard limits : ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(~isfield(r.softlims.RATE_A, 'overload'), [t 'no softlims.RATE_A.overload']);
t_ok(~isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'no softlims.RATE_A.ovl_cost']);
t_is(r.f, 9513.613051, 4, [t 'f']);
t_is(r.gen(:, PG), [17.759246; 65.641269; 240], 4, [t 'Pg']);
t_is(r.branch(:, PF), [17.759246; -25.216499; -115.346644; 240; 0; 120; 18.577311; -65.641269; 84.170657; -42.784248], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-17.759246; 25.346644; 120; -240; 0; -118.577311; -18.529387; 65.641269; -82.215752; 42.975745], 4, [t 'Pt']);
t_is(r.branch(:, QF), [15.646152; 6.577259; -6.087866; 8.021258; 0; -4.520662; -26.623563; -6.687069; -2.629290; -27.085546], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-15.370220; -23.912134; -15.548682; 20.069344; 0; -8.376437; 9.316359; 8.954090; -22.914454; 8.792961], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 29.248033; 0; 0; 8.417115; 0; 0; 0; 0], 4, [t 'mu Sf+St']);

t = 'AC (flow_lim ''P'') - soft limits (satisfied) : ';
mpc = toggle_softlims(mpc, 'on');
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9513.613051, 4, [t 'f']);
t_is(r.gen(:, PG), [17.759246; 65.641269; 240], 4, [t 'Pg']);
t_is(r.branch(:, PF), [17.759246; -25.216499; -115.346644; 240; 0; 120; 18.577311; -65.641269; 84.170657; -42.784248], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-17.759246; 25.346644; 120; -240; 0; -118.577311; -18.529387; 65.641269; -82.215752; 42.975745], 4, [t 'Pt']);
t_is(r.branch(:, QF), [15.646152; 6.577259; -6.087866; 8.021258; 0; -4.520662; -26.623563; -6.687069; -2.629290; -27.085546], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-15.370220; -23.912134; -15.548682; 20.069344; 0; -8.376437; 9.316359; 8.954090; -22.914454; 8.792961], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 29.248033; 0; 0; 8.417115; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.RATE_A.overload, zeros(10, 1), 12, [t 'softlims.RATE_A.overload']);
t_is(r.softlims.RATE_A.ovl_cost, zeros(10, 1), 12, [t 'softlims.RATE_A.ovl_cost']);
t_is(r.order.branch.status.on(r.order.int.softlims.RATE_A.idx), [3; 6; 8; 9; 10], 12, [t 'softlims.RATE_A.idx']);

t = 'AC (flow_lim ''P'') - soft limits (violated) : ';
mpc = rmfield(mpc, 'softlims');
mpc.softlims.RATE_A.type = 'unbnd';
mpc.softlims.RATE_A.cost = 20;
for prop = {'VMIN', 'VMAX', 'ANGMIN', 'ANGMAX', 'PMIN', 'PMAX', 'QMIN', 'QMAX'}
    prop = prop{1};
    mpc.softlims.(prop).type = 'none';
end
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9476.854715, 4, [t 'f']);
t_is(r.gen(:, PG), [10; 69.833910; 244.077740; ], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -28.933383; -119.111239; 244.077740; 0; 120; 18.578292; -69.833910; 88.362736; -38.762267], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-10; 29.111239; 124.077740; -244.077740; 0; -118.578292; -18.528826; 69.833910; -86.237733; 38.933383], 4, [t 'Pt']);
t_is(r.branch(:, QF), [22.204979; 10.325903; -2.434897; 5.980065; 0; -5.159694; -27.240912; -4.815545; -5.105567; -30.250169], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-21.917920; -27.565103; -17.970864; 23.130558; 0; -7.759088; 9.921112; 7.362539; -19.749831; 11.592017], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 20; 0; 0; 11.526783; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.RATE_A.overload, [0; 0; 4.07774; 0; 0; 0; 0; 0; 0; 0], 6, [t 'softlims.RATE_A.overload']);
t_is(r.softlims.RATE_A.ovl_cost, [0; 0; 81.554809; 0; 0; 0; 0; 0; 0; 0], 4, [t 'softlims.RATE_A.ovl_cost']);



if have_fcn('matlab')
    warning(s2.state, sing_matrix_warn_id);
end

t_end;
%%
function overload_loop(r, t, yn, lst)
%%% r is the result structure
%%% t is the test string
%%% yn is whether to check if it is (1) or is not (0) a field
%%% lst is a cell of properties
%%%     if yn == 1 i.e. checking if properties exist then this is a list of
%%%     properties to check
%%%     if yn == 0 i.e. checking that the properties do not exist, this is
%%%     a list of properties to exclude

if nargin == 3
    lst = {};
end

props_default = {'RATE_A', 'VMAX', 'VMIN', 'PMAX', 'PMIN', 'QMAX', 'QMIN', 'ANGMAX', 'ANGMIN'};

if ~yn
    lst = props_default(~ismember(props_default, lst));
else
    if isempty(lst)
        lst = props_default;
    else
        lst = props_default(ismember(props_default, lst));
    end
end
for p = lst
    p = p{1};
    if yn
        t_ok(isfield(r.softlims.(p), 'overload'), [t sprintf('softlims.%s.overload exists',p)]);
        t_ok(isfield(r.softlims.(p), 'ovl_cost'), [t sprintf('softlims.%s.ovl_cost exists',p)]);
    else
        t_ok(~isfield(r.softlims.(p), 'overload'), [t sprintf('no softlims.%s.overload',p)]);
        t_ok(~isfield(r.softlims.(p), 'ovl_cost'), [t sprintf('no softlims.%s.ovl_cost',p)]);
    end
end

function r = toggle_run_check(mpc, mpopt, t, on_off)
%%% checks the softlims status of mpc, runs the opf, checks for success and
%%% returns the results structure

if ~on_off
    t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
else
    t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
end
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);

function mu_cost_test(r,t)
%%% Tests whether the shadow price on violations matches the respective
%%% softlim cost. The shadow price should equal the softlim cost when the
%%% slack variable is used but the constraint is not binding. i.e. there is
%%% a non-zero overload but the slack variable is not at its own limit.
%%% the function therefore searches for any non-zero overload (mumask1)
%%% that is also below the slack variables upper bound stored in s.ub
%%% (mumask2). The cost of these variables is compared to the respective
%%% shadow prices.

[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
lims = struct('VMAX', 'bus', 'VMIN','bus', 'ANGMAX', 'branch', 'ANGMIN', 'branch',...
							'RATE_A', 'branch', 'PMAX', 'gen', 'PMIN', 'gen', 'QMAX', 'gen', 'QMIN', 'gen'); 
for prop = fieldnames(lims).'
    prop = prop{1};
    mat  = lims.(prop);
    if ~isfield(r.softlims,prop)
        continue
    end
    if strcmp(r.softlims.(prop).type, 'none')
        continue
    end
    s = r.softlims.(prop);
    % slack variable is non zero
    mumask1 = s.overload > 1e-6; 
    % ensure that overloads are not at the slack variable boundary since
    % then the shadow price can take on any value again.
    mumask2 = true(size(mumask1));
    mumask2(s.idx) = s.overload(s.idx) < s.ub;
    
    mumask = mumask1 & mumask2;
    
    % linear cost
    cst = zeros(size(mumask1));
    cst(s.idx) = s.cost;
    cst = cst(mumask); %keep only relevant entries
    
    if strcmp(prop,'RATE_A')
        mu = sum(r.branch(mumask,MU_SF:MU_ST));
    else
        muidx = eval(['MU_' prop]);
        mu = r.(mat)(mumask,muidx);       
    end
    
    t_is(mu,cst,4,[t 'mu ' prop '=cost'])
end