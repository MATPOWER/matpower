function t_opf_softlims(quiet)
%T_OPF_SOFTLIMS  Tests for userfcn callbacks (softlims) w/OPF.
%   Includes high-level tests of soft limits implementations.

%   MATPOWER
%   Copyright (c) 2009-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Eran Schweitzer, Arizona State University
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

casefile = 'case9';
fname_ac = 'pretty_print_softlims_ac.txt';
fname_dc = 'pretty_print_softlims_dc.txt';
fname = 'pretty_print_softlims';
rn = fix(1e9*rand);
tmp_fname_ac = sprintf('%s_ac_%d.txt', fname, rn);
tmp_fname_dc = sprintf('%s_dc_%d.txt', fname, rn);

if quiet
    verbose = 0;
else
    verbose = 0;
end

t_begin(842, quiet);

%% define constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% set options
%mpopt = mpoption('opf.dc.solver', 'MIPS', 'opf.ac.solver', 'IPOPT');
mpopt = mpoption('opf.dc.solver', 'MIPS', 'opf.ac.solver', 'MIPS');
%mpopt = mpoption('opf.dc.solver', 'GLPK', 'opf.ac.solver', 'FMINCON');
%mpopt = mpoption('opf.dc.solver', 'GUROBI', 'opf.ac.solver', 'IPOPT');
%mpopt = mpoption('opf.dc.solver', 'CPLEX', 'opf.ac.solver', 'KNITRO');
%mpopt = mpoption(mpopt, 'fmincon.tol_x', 1e-9, 'fmincon.tol_f', 1e-9);
%mpopt = mpoption(mpopt, 'ipopt.opts.tol', 1e-10, 'ipopt.opts.acceptable_tol', 1e-10);
%mpopt = mpoption(mpopt, 'knitro.tol_x', 1e-9, 'knitro.tol_f', 1e-9);
mpopt = mpoption(mpopt, 'opf.violation', 1e-6, 'mips.gradtol', 1e-8, ...
        'mips.comptol', 1e-8, 'mips.costtol', 1e-9);
if verbose <= 1
    mpopt = mpoption(mpopt, 'out.all', 0);
end
mpopt = mpoption(mpopt, 'verbose', verbose);

if have_fcn('octave')
    sing_matrix_warn_id = 'Octave:singular-matrix';
    sing_matrix_warn_id = 'Octave:nearly-singular-matrix';
    file_in_path_warn_id = 'Octave:data-file-in-path';
    s1 = warning('query', file_in_path_warn_id);
    warning('off', file_in_path_warn_id);
else
    sing_matrix_warn_id = 'MATLAB:singularMatrix';
    sing_matrix_warn_id = 'MATLAB:nearlySingularMatrix';
end
s2 = warning('query', sing_matrix_warn_id);
warning('off', sing_matrix_warn_id);

%% load and modify case file
mpc = loadcase(casefile);
mpc.bus(:,BUS_I) = 10*mpc.bus(:,BUS_I);
mpc.branch(:,[F_BUS,T_BUS]) = 10*mpc.branch(:,[F_BUS,T_BUS]);
mpc.gen(:,GEN_BUS) = 10*mpc.gen(:,GEN_BUS);
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

%% reorder generators so they are not consecutive
mpc.gen = mpc.gen([3,1,2],:);
mpc.gencost = mpc.gencost([3,1,2],:);
%% duplicate 3rd generator and make it off-line
mpc.gen = [mpc.gen(1,:); mpc.gen(3,:); mpc.gen(2:end,:)];
mpc.gencost = [mpc.gencost(1,:); mpc.gencost(3,:); mpc.gencost(2:end,:)];
mpc.gen(2,GEN_STATUS) = 0;
%% create soft limit inputs (emulates previous implementation just for branches)
mpc.softlims.RATE_A.hl_mod = 'remove';
mpc.softlims.RATE_A.idx  = (2:nl)';
mpc.softlims.RATE_A.cost = 100 * ones(nl-1, 1);
for lim = {'VMIN', 'VMAX', 'ANGMIN', 'ANGMAX', 'PMIN', 'PMAX', 'QMIN', 'QMAX'}
    mpc.softlims.(lim{:}).hl_mod = 'none';
end
mpc0 = mpc;     %% save to initialize later runs

%% mixed softlims structure: uses different types of limits
sdefault = struct();
for lm = {'RATE_A', 'VMIN', 'VMAX', 'ANGMIN', 'ANGMAX', 'PMIN', 'PMAX', ...
        'QMIN', 'QMAX'}
    lim = lm{:};
    switch lim
        case 'RATE_A'
            sdefault.(lim).hl_mod = 'scale';
            sdefault.(lim).hl_val = 1.5;
        case 'VMIN'
            sdefault.(lim).hl_mod = 'scale';
            sdefault.(lim).hl_val = 0.5;
        case 'VMAX'
            sdefault.(lim).hl_mod = 'scale';
            sdefault.(lim).hl_val = 1.5;
        case {'PMAX', 'QMAX', 'QMIN'}
            sdefault.(lim).hl_mod = 'remove';
        case 'PMIN'
            sdefault.(lim).hl_mod = 'scale';
            sdefault.(lim).hl_val = 0;
        case 'ANGMAX'
            sdefault.(lim).hl_mod = 'replace';
            sdefault.(lim).hl_val = 360;
        case 'ANGMIN'
            sdefault.(lim).hl_mod = 'replace';
            sdefault.(lim).hl_val = -360;
    end
end
%% generator ordering
mpc = mpc0;
mpc.bus(:,QD) = mpc.bus(:,QD)*2;

mpc.gen(:,[QMAX,QMIN]) = [1.1,-1.1; 2.2,-2.2; 3.3,-3.3; 4.4, -4.4];
mpc.gen(:,PMIN) = [1.1; 2.2; 3.3; 4.4];

% mpc.softlims = sdefault;
mpc.softlims.QMAX = rmfield(mpc.softlims.QMAX, 'hl_mod');
mpc.softlims.QMAX.idx = [1,3];

t = 'generator ordering: ';
mpc = toggle_softlims(mpc,'on');
r = toggle_run_check(mpc, mpopt, t, 1);
gen_order_check(mpc, r, t, mpopt);
% opf.softlims.default = 0
mpopt.opf.softlims.default = 0;
mpc = mpc0;
mpc = rmfield(mpc,'softlims');

t = 'mpopt - opf.softlims.default = 0 (all none): ';
mpc = toggle_softlims(mpc,'on');
r = toggle_run_check(mpc, mpopt, t, 1);
for lm = fieldnames(r.softlims).'
    lim = lm{:};
    t_ok(strcmp(r.softlims.(lim).hl_mod, 'none'), [t lim '.hl_mod = ''none'''])
end
gen_order_check(mpc, r, t, mpopt);

t = 'mpopt - opf.softlims.default = 0 (w/RATE_A): ';
mpc.softlims.RATE_A = struct('hl_mod', 'scale', 'hl_val', 1.5); %initializing an empty strucure results in default
r = toggle_run_check(mpc, mpopt, t, 1);
vv = r.om.get_idx();
[x0, xmin, xmax] = r.om.params_var();
for lm = fieldnames(r.softlims).'
    lim = lm{:};
    if ~strcmp(lim, 'RATE_A')
        t_ok(strcmp(r.softlims.(lim).hl_mod, 'none'), [t lim '.hl_mod = ''none'''])
    else
        t_ok(strcmp(r.softlims.RATE_A.hl_mod,'scale'), [t lim '.hl_mod = ''scale'''])
        varname = ['s_' lower(lim)];
        ub = xmax(vv.i1.(varname):vv.iN.(varname)) * r.baseMVA;
        t_is(ub, 0.5*r.branch(r.softlims.RATE_A.idx,RATE_A), 5, [t '.hl_val = 0.5'])
    end
end
gen_order_check(mpc, r, t, mpopt);

mpopt.opf.softlims.default = 1;
mpc = mpc0;
mpc = rmfield(mpc,'softlims');
% mpc.softlims = sdefault;
t = 'AC - opf.softlims.default = 1 (all default): ';
mpc = toggle_softlims(mpc,'on');
r = toggle_run_check(mpc, mpopt, t, 1);
vv = r.om.get_idx();
[x0, xmin, xmax] = r.om.params_var();
for lm = fieldnames(r.softlims).'
    lim = lm{:};
    varname = ['s_' lower(lim)];
    ub = xmax(vv.i1.(varname):vv.iN.(varname));
    switch lim
        case {'RATE_A', 'PMAX', 'PMIN', 'QMAX', 'QMIN'}
            ub = ub * r.baseMVA;
        case {'ANGMAX', 'ANGMIN'}
            ub = ub * 180/pi;
    end
    switch lim
        case {'PMAX', 'QMAX', 'RATE_A', 'ANGMAX', 'VMAX'}
            %% default hl_mod = 'remove'
            t_ok(strcmp(r.softlims.(lim).hl_mod,'remove'), [t lim '.hl_mod = ''remove'''])
            t_ok(r.softlims.(lim).hl_val == Inf, [t lim '.hl_val = Inf'])
            t_ok(all(isinf(ub)), [t lim ' ub = Inf'])
        case {'QMIN', 'ANGMIN'}
            %% default hl_mod = 'remove'
            t_ok(strcmp(r.softlims.(lim).hl_mod,'remove'), [t lim '.hl_mod = ''remove'''])
            t_ok(r.softlims.(lim).hl_val == -Inf, [t lim '.hl_val = -Inf'])
            t_ok(all(isinf(ub)), [t lim ' ub = Inf'])
        case 'VMIN'
            %% default hl_mod = 'replace'
            t_ok(strcmp(r.softlims.(lim).hl_mod,'replace'), [t lim '.hl_mod = ''replace'''])
            t_ok(r.softlims.(lim).hl_val == 0, [t lim '.hl_val = 0'])
            t_ok(all(ub == r.bus(:,VMIN)), [t lim ' ub = VMIN'])
        case 'PMIN'
            %% default hl_mod = 'replace'
            t_ok(strcmp(r.softlims.(lim).hl_mod,'replace'), [t lim '.hl_mod = ''replace'''])
            t_ok(all(r.softlims.(lim).hl_val(r.gen(r.softlims.PMIN.idx, PMIN) > 0) == 0), [t lim '.hl_val = 0 (gens)'])
            t_ok(all(ub(r.gen(r.softlims.PMIN.idx, PMIN) > 0) == r.gen(r.softlims.PMIN.idx,PMIN)), [t lim ' ub=PMIN (gens)'])
    end
end
gen_order_check(mpc, r, t, mpopt);

t = 'DC - opf.softlims.default = 1 (all default): ';
mpopt1 = mpoption(mpopt, 'model', 'DC');
r = toggle_run_check(mpc, mpopt1, t, 1);
vv = r.om.get_idx();
[x0, xmin, xmax] = r.om.params_var();
for lm = fieldnames(r.softlims).'
    lim = lm{:};
    varname = ['s_' lower(lim)];
    ub = xmax(vv.i1.(varname):vv.iN.(varname));
    switch lim
        case {'RATE_A', 'PMAX', 'PMIN', 'QMAX', 'QMIN'}
            ub = ub * r.baseMVA;
        case {'ANGMAX', 'ANGMIN'}
            ub = ub * 180/pi;
    end
    switch lim
        case {'PMAX', 'QMAX', 'RATE_A', 'ANGMAX', 'VMAX'}
            %% default hl_mod = 'remove'
            t_ok(strcmp(r.softlims.(lim).hl_mod,'remove'), [t lim '.hl_mod = ''remove'''])
            t_ok(r.softlims.(lim).hl_val == Inf, [t lim '.hl_val = Inf'])
            t_ok(all(isinf(ub)), [t lim ' ub = Inf'])
        case {'QMIN', 'ANGMIN'}
            %% default hl_mod = 'remove'
            t_ok(strcmp(r.softlims.(lim).hl_mod,'remove'), [t lim '.hl_mod = ''remove'''])
            t_ok(r.softlims.(lim).hl_val == -Inf, [t lim '.hl_val = -Inf'])
            t_ok(all(isinf(ub)), [t lim ' ub = Inf'])
        case 'VMIN'
            %% default hl_mod = 'replace'
            t_ok(strcmp(r.softlims.(lim).hl_mod,'replace'), [t lim '.hl_mod = ''replace'''])
            t_ok(r.softlims.(lim).hl_val == 0, [t lim '.hl_val = 0'])
            t_ok(all(ub == r.bus(:,VMIN)), [t lim ' ub = VMIN'])
        case 'PMIN'
            %% default hl_mod = 'replace'
            t_ok(strcmp(r.softlims.(lim).hl_mod,'replace'), [t lim '.hl_mod = ''replace'''])
            t_ok(all(r.softlims.(lim).hl_val(r.gen(r.softlims.PMIN.idx, PMIN) > 0) == 0), [t lim '.hl_val = 0 (gens)'])
            t_ok(all(ub(r.gen(r.softlims.PMIN.idx, PMIN) > 0) == r.gen(r.softlims.PMIN.idx,PMIN)), [t lim ' ub=PMIN (gens)'])
    end
end
gen_order_check(mpc, r, t, mpopt1);

%% Default settings check
mpc = mpc0;
mpc.softlims = sdefault;

t = 'Defaults - no active softlims: ';
r = toggle_run_check(mpc, mpopt, t, 0);
delta = calc_branch_angle(r);
overload_loop(r, t, 0);
t_is(r.bus(:,VM), [1.040963;1.1;1.1;1.042513;1.041787;1.096141;1.079794;1.089516;1.031234], 4, [t 'Vm']);
t_is(r.branch(:,PF), [17.3121  -25.0634 -115.1668  239.9012         0  119.9646   18.5301  -66.4562   84.9411  -42.2020].', 4, [t 'PF'])
t_is(r.branch(:,QF), [-2.7212   -2.8421  -16.2418   21.2760         0   -2.9141  -25.3254  -17.1234    8.9449  -17.4049].', 4, [t 'QF'])
t_is(r.branch(:,PT), [-17.3121   25.1668  119.9366 -239.9012         0 -118.5301  -18.4848   66.4562  -82.7980   42.3756].',4, [t 'PT'])
t_is(r.branch(:,QT), [2.8844  -13.7582   -3.9017    6.8158         0   -9.6746    8.1785   19.6031  -32.5951   -0.0423].', 4, [t 'QT'])
t_is(delta, [0.5265   -1.2681   -9.9352    6.6955         0    5.8081    0.7187   -1.9861    6.5458   -1.8692].', 4, [t 'delta'])
t_is(r.gen(:,PG), [239.9012         0   17.3121   66.4562].', 4, [t 'PG'])
t_is(r.gen(:,QG), [21.2760         0   -2.7212   19.6031].', 4, [t 'QG'])
t_is(r.bus(:,MU_VMAX), [0;526.455671;232.322240;0;0;0;0;0;0], 4, [t 'mu VM ub'])
t_is(r.bus(:,MU_VMIN), zeros(9,1), 4, [t 'mu VM lb'])
t_is(sum(r.branch(:,MU_SF:MU_ST),2), [0.0000         0   29.2428         0         0    8.5610         0    0.0000    0.0000    0.0000].', 4, [t 'mu SF+ST'])
t_is(r.branch(:,MU_ANGMAX), zeros(10,1), 4, [t 'mu ANGMAX'])
t_is(r.branch(:,MU_ANGMIN), zeros(10,1), 4, [t 'mu ANGMIN'])
t_is(r.gen(:,MU_PMAX), zeros(4,1), 4, [t 'mu PMAX'])
t_is(r.gen(:,MU_PMIN), zeros(4,1), 4, [t 'mu PMIN'])
t_is(r.gen(:,MU_QMAX), zeros(4,1), 4, [t 'mu QMAX'])
t_is(r.gen(:,MU_QMIN), zeros(4,1), 4, [t 'mu QMIN'])
gen_order_check(mpc, r, t, mpopt);

t = 'Defaults - softlims on: ';
mpc = toggle_softlims(mpc,'on');
r = toggle_run_check(mpc, mpopt, t, 1);
delta = calc_branch_angle(r);
overload_loop(r, t, 1);
t_is(r.bus(:,VM), [1.040963;1.1;1.1;1.042513;1.041787;1.096141;1.079794;1.089516;1.031234], 4, [t 'Vm']);
t_is(r.branch(:,PF), [17.3121  -25.0634 -115.1668  239.9012         0  119.9646   18.5301  -66.4562   84.9411  -42.2020].', 4, [t 'PF'])
t_is(r.branch(:,QF), [-2.7212   -2.8421  -16.2418   21.2760         0   -2.9141  -25.3254  -17.1234    8.9449  -17.4049].', 4, [t 'QF'])
t_is(r.branch(:,PT), [-17.3121   25.1668  119.9366 -239.9012         0 -118.5301  -18.4848   66.4562  -82.7980   42.3756].',4, [t 'PT'])
t_is(r.branch(:,QT), [2.8844  -13.7582   -3.9017    6.8158         0   -9.6746    8.1785   19.6031  -32.5951   -0.0423].', 4, [t 'QT'])
t_is(delta, [0.5265   -1.2681   -9.9352    6.6955         0    5.8081    0.7187   -1.9861    6.5458   -1.8692].', 4, [t 'delta'])
t_is(r.gen(:,PG), [239.9012         0   17.3121   66.4562].', 4, [t 'PG'])
t_is(r.gen(:,QG), [21.2760         0   -2.7212   19.6031].', 4, [t 'QG'])
t_is(r.bus(:,MU_VMAX), [0;526.455671;232.322240;0;0;0;0;0;0], 4, [t 'mu VM ub'])
t_is(r.bus(:,MU_VMIN), zeros(9,1), 4, [t 'mu VM lb'])
t_is(sum(r.branch(:,MU_SF:MU_ST),2), [0.0000         0   29.2428         0         0    8.5610         0    0.0000    0.0000    0.0000].', 4, [t 'mu SF+ST'])
t_is(r.branch(:,MU_ANGMAX), zeros(10,1), 4, [t 'mu ANGMAX'])
t_is(r.branch(:,MU_ANGMIN), zeros(10,1), 4, [t 'mu ANGMIN'])
t_is(r.gen(:,MU_PMAX), zeros(4,1), 4, [t 'mu PMAX'])
t_is(r.gen(:,MU_PMIN), zeros(4,1), 4, [t 'mu PMIN'])
t_is(r.gen(:,MU_QMAX), zeros(4,1), 4, [t 'mu QMAX'])
t_is(r.gen(:,MU_QMIN), zeros(4,1), 4, [t 'mu QMIN'])
gen_order_check(mpc, r, t, mpopt);

%% tighten limits to get violations
t = 'Defaults - softlims w/overloads: ';
mpc = mpc0;
mpc.softlims = sdefault;
mpc.bus(4,VMAX) = 0.9;
mpc.softlims.VMAX.cost = 75000;
mpc.bus(9,VMIN)  = 1.05;
mpc.softlims.VMIN.cost = 75000;
mpc.softlims.VMIN.busnum = 10*[1 2 8 9]';
mpc.softlims.RATE_A.cost = 10;
mpc.gen(3,PMIN) = 25; 
mpc.softlims.PMIN.cost = 5;
mpc.gen(1,PMAX) = 225; 
mpc.softlims.PMAX.idx  = [1;2;3;4];
mpc.softlims.PMAX.cost = [3;1000;1000;1000];
mpc.gen(4,QMAX) = 50; 
mpc.softlims.QMAX.idx  = [1;2;3;4];
mpc.softlims.QMAX.cost = [1000;1000;1000;20];
mpc.gen(3,QMIN) = -45; 
mpc.softlims.QMIN.idx  = [1;2;3;4];
mpc.softlims.QMIN.cost = [1000; 1000; 10; 1000];
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
mu_cost_test(r, t, mpopt);
gen_order_check(mpc, r, t, mpopt);

t = 'printpf - AC';
[fd, msg] = fopen(tmp_fname_ac, 'at');
if fd == -1
    error(msg);
else
    mpopt2 = mpoption(mpopt, 'out.all', -1, 'out.bus', 0, 'out.branch', 0, 'out.sys_sum', 0);
    printpf(r, fd, mpopt2);
    fclose(fd);
end
got = fileread(tmp_fname_ac);
got = strrep(got, char([13 10]), char(10));             %% Win to Unix EOL chars
got = regexprep(got, 'Converged in (.*) seconds', 'Converged in 0.00 seconds');
expected = fileread(fname_ac);
expected = strrep(expected, char([13 10]), char(10));   %% Win to Unix EOL chars
t_ok(strcmp(got, expected), t);
delete(tmp_fname_ac);

t = 'Repeat solve: ';
r2 = runopf(r, mpopt);
t_is(r2.bus, r.bus, 4, [t 'bus matrices are the same'])
t_is(r2.branch, r.branch, 4, [t 'branch matrices are the same'])
t_is(r2.gen, r.gen, 4, [t 'bus matrices are the same'])

%% test save case
t = 'savecase(fname, r) : ';
fn = sprintf('softlims_savecase_test_%d', fix(1e8*rand));
savecase(fn, r);
mpc1 = loadcase(fn);
delete([fn '.m']);
t_ok(isfield(mpc1, 'softlims'), [t 'mpc.softlims']);
for lm = fieldnames(r.softlims).'
    lim = lm{:};
    t_ok(isfield(mpc1.softlims, lim), [t 'mpc.softlims.' lim])
    for f = fieldnames(r.softlims.(lim)).'
        field = f{:};
        if any(isinf(r.softlims.(lim).(field)))
            t_ok(isequal(mpc1.softlims.(lim).(field), r.softlims.(lim).(field)), [t 'mpc.softlims.' lim '.' field])
        else
            t_is(mpc1.softlims.(lim).(field), r.softlims.(lim).(field), 5, [t 'mpc.softlims.' lim '.' field])
        end
    end
end

%% defaults with DC
mpc = mpc0;
mpc.softlims = sdefault;
t = 'DC - hard limits : ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
mpopt1 = mpoption(mpopt, 'model', 'DC');
r = rundcopf(mpc, mpopt1);
t_ok(r.success, [t 'success']);
delta = calc_branch_angle(r);
overload_loop(r, t, 0, {}, 1);
t_is(r.gen(:, PG), [240.0000         0   12.6870   62.3130].', 4, [t 'Pg']);
t_is(r.branch(:, PF), [12.687; -30; -120; 240; 0; 120; 20; -62.3130; 82.3130; -42.687], 4, [t 'Pf']);
t_is(delta, [0.4187   -1.5814  -11.6883    8.0581         0    6.9305    0.8251   -2.2314    7.5931   -2.0789].', 4, [t 'delta'])
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 35.6504; 0; 0; 7.9756; 0; 0; 0; 0], 4, [t 'mu Pf']);
gen_order_check(mpc, r, t, mpopt1);

t = 'DC - softlims on: ';
mpc = toggle_softlims(mpc,'on');
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = rundcopf(mpc, mpopt1);
t_ok(r.success, [t 'success']);
delta = calc_branch_angle(r);
overload_loop(r, t, 1, {}, 1);
t_is(r.gen(:, PG), [240.0000         0   12.6870   62.3130].', 4, [t 'Pg']);
t_is(r.branch(:, PF), [12.687; -30; -120; 240; 0; 120; 20; -62.3130; 82.3130; -42.687], 4, [t 'Pf']);
t_is(delta, [0.4187   -1.5814  -11.6883    8.0581         0    6.9305    0.8251   -2.2314    7.5931   -2.0789].', 4, [t 'delta'])
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 35.6504; 0; 0; 7.9756; 0; 0; 0; 0], 4, [t 'mu Pf']);
gen_order_check(mpc, r, t, mpopt1);

t = 'DC - regression : ';
mpc = mpc0;
mpc = toggle_softlims(mpc,'on');
t_ok(toggle_softlims(mpc, 'status'), [t 'toggle_softlims(mpc, ''status'') == 1']);
mpc.softlims = struct();
mpc.softlims.RATE_A.hl_mod = 'none';
r = rundcopf(mpc, mpopt1);
t_ok(r.success, [t 'success']);

%% tighten limits to get violations
t = 'DC - softlims with overloads: ';
mpc = mpc0;
mpc.softlims = sdefault;
mpc.branch([3,6,7,8],RATE_A) = 75;
mpc.softlims.RATE_A.cost = 15;
mpc.gen(3,PMIN) = 75; 
mpc.gen(4,PMIN) = 80;
mpc.softlims.PMIN.idx = [1;2;3;4];
mpc.softlims.PMIN.cost = [1000;1000;5;4];
mpc.gen(1,PMAX) = 90; 
mpc.softlims.PMAX.idx = [1;2;3;4];
mpc.softlims.PMAX.cost = [10;1000;1000;1000];
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
r = rundcopf(mpc, mpopt1);
t_ok(r.success, [t 'success']);
overload_loop(r, t, 1, {}, 1);
mu_cost_test(r, t, mpopt1);
gen_order_check(mpc, r, t, mpopt1);

t = 'printpf - DC';
[fd, msg] = fopen(tmp_fname_dc, 'at');
if fd == -1
    error(msg);
else
    mpopt2 = mpoption(mpopt1, 'out.all', -1, 'out.bus', 0, 'out.branch', 0, 'out.sys_sum', 0);
    printpf(r, fd, mpopt2);
    fclose(fd);
end
got = fileread(tmp_fname_dc);
got = strrep(got, char([13 10]), char(10));             %% Win to Unix EOL chars
got = regexprep(got, 'Converged in (.*) seconds', 'Converged in 0.00 seconds');
expected = fileread(fname_dc);
expected = strrep(expected, char([13 10]), char(10));   %% Win to Unix EOL chars
t_ok(strcmp(got, expected), t);
delete(tmp_fname_dc);

%% unbounded limits
mpc = mpc0;
mpc.softlims.RATE_A.hl_mod = 'none';
mpc.softlims.VMAX.hl_mod   = 'remove';
mpc.softlims.VMIN.hl_mod   = 'remove';
t = 'V lims - hard limits: ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
overload_loop(r, t, 0);
t_is(r.bus(:,VM), [1.040963;1.1;1.1;1.042513;1.041787;1.096141;1.079794;1.089516;1.031234], 4, [t 'Vm']);
t_is(r.bus(:,MU_VMAX), [0;526.455671;232.322240;0;0;0;0;0;0], 4, [t 'mu VM ub'])
gen_order_check(mpc, r, t, mpopt);

mpc = toggle_softlims(mpc,'on');
t = 'V lims (remove) - softlims (satisfied): ';
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.bus(:,VM), [1.040963;1.1;1.1;1.042513;1.041787;1.096141;1.079794;1.089516;1.031234], 4, [t 'Vm']);
t_is(r.bus(:,MU_VMAX), [0;526.455671;232.322240;0;0;0;0;0;0], 4, [t 'mu VM ub'])
gen_order_check(mpc, r, t, mpopt);

mpc.bus(4:9, [VMIN, VMAX]) = 1;
t = 'V lims (remove) - softlims w/overloads: ';
r = toggle_run_check(mpc, mpopt, t, 1);
overload_loop(r, t, 1, {'VMIN', 'VMAX'})
overload_loop(r, t, 0, {'VMIN', 'VMAX'})
mu_cost_test(r, t, mpopt);
gen_order_check(mpc, r, t, mpopt);

%% scale modification
mpc = mpc0;
mpc.softlims = sdefault;
for lm = fieldnames(mpc.softlims).'
    lim = lm{:};
    if ~ismember(lim,{'VMIN','VMAX'})
        mpc.softlims.(lim).hl_mod = 'none';
    end
end
mpc.bus(4:9, [VMIN, VMAX]) = 1;
mpc = toggle_softlims(mpc,'on');
t = 'V lims (scale) - softlims w/overloads: ';
r = toggle_run_check(mpc, mpopt, t, 1);
overload_loop(r, t, 1, {'VMIN', 'VMAX'})
overload_loop(r, t, 0, {'VMIN', 'VMAX'})
mu_cost_test(r, t, mpopt);
gen_order_check(mpc, r, t, mpopt);

%% replace modification
mpc = mpc0;
mpc.softlims = sdefault;
for lm = fieldnames(mpc.softlims).'
    lim = lm{:};
    if strcmp(lim,'VMIN')
        mpc.softlims.VMIN.hl_mod = 'replace';
        mpc.softlims.VMIN.hl_val   = 0.7;
    elseif strcmp(lim,'VMAX')
        mpc.softlims.VMAX.hl_mod = 'replace';
        mpc.softlims.VMAX.hl_val   = 1.2;
    else
        mpc.softlims.(lim).hl_mod = 'none';
    end
end
mpc.bus(4:9, [VMIN, VMAX]) = 1;
mpc = toggle_softlims(mpc,'on');
t = 'V lims (replace) - softlims w/overloads: ';
r = toggle_run_check(mpc, mpopt, t, 1);
overload_loop(r, t, 1, {'VMIN', 'VMAX'})
overload_loop(r, t, 0, {'VMIN', 'VMAX'})
vv = r.om.get_idx();
[x0, xmin, xmax] = r.om.params_var();
ub = xmax(vv.i1.s_vmin:vv.iN.s_vmin);
t_ok(all(ub == (r.bus(r.softlims.VMIN.idx,VMIN) - 0.7)), [t 'VMIN ub = VMIN - 0.7'])
ub = xmax(vv.i1.s_vmax:vv.iN.s_vmax);
t_ok(all(ub == (1.2 - r.bus(r.softlims.VMAX.idx,VMAX))), [t 'VMIN ub = 1.2 - VMAX'])
mu_cost_test(r, t, mpopt);
gen_order_check(mpc, r, t, mpopt);

%% shift modification
mpc = mpc0;
mpc.softlims = sdefault;
for lm = fieldnames(mpc.softlims).'
    lim = lm{:};
    if strcmp(lim,'VMIN')
        mpc.softlims.VMIN.hl_mod = 'shift';
        mpc.softlims.VMIN.hl_val   = 0.2;
    elseif strcmp(lim,'VMAX')
        mpc.softlims.VMAX.hl_mod = 'shift';
        mpc.softlims.VMAX.hl_val   = 0.2;
    else
        mpc.softlims.(lim).hl_mod = 'none';
    end
end
mpc.bus(4:9, [VMIN, VMAX]) = 1;
mpc = toggle_softlims(mpc,'on');
t = 'V lims (shift) - softlims w/overloads: ';
r = toggle_run_check(mpc, mpopt, t, 1);
overload_loop(r, t, 1, {'VMIN', 'VMAX'})
overload_loop(r, t, 0, {'VMIN', 'VMAX'})
vv = r.om.get_idx();
[x0, xmin, xmax] = r.om.params_var();
ub = xmax(vv.i1.s_vmin:vv.iN.s_vmin);
t_ok(all(ub == 0.2), [t 'VMIN ub = shift (0.2)'])
ub = xmax(vv.i1.s_vmax:vv.iN.s_vmax);
t_ok(all(ub == 0.2), [t 'VMAX ub = shift (0.2)'])
mu_cost_test(r, t, mpopt);
gen_order_check(mpc, r, t, mpopt);
%% Tests
%% the following are all the test from the original softlims implementation
%% that only included RATE_A. The only changes are adaptations to the
%% changed data structure.
mpc = mpc0;
t = 'DC - hard limits : ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
mpopt1 = mpoption(mpopt, 'model', 'DC');
r = rundcopf(mpc, mpopt1);
t_ok(r.success, [t 'success']);
t_ok(~isfield(r.softlims.RATE_A, 'overload'), [t 'no softlims.RATE_A.overload']);
t_ok(~isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'no softlims.RATE_A.ovl_cost']);
t_is(r.f, 9126.87, 4, [t 'f']);
t_is(r.gen(:, PG), [240.0000         0   12.6870   62.3130].', 4, [t 'Pg']);
t_is(r.branch(:, PF), [12.687; -30; -120; 240; 0; 120; 20; -62.3130; 82.3130; -42.687], 4, [t 'Pf']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 35.6504; 0; 0; 7.9756; 0; 0; 0; 0], 4, [t 'mu Pf']);
gen_order_check(mpc, r, t, mpopt1);

t = 'DC - softlims (satisfied) : ';
mpc = toggle_softlims(mpc, 'on');
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = rundcopf(mpc, mpopt1);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9126.87, 4, [t 'f']);
t_is(r.gen(:, PG), [240.0000         0   12.6870   62.3130].', 4, [t 'Pg']);
t_is(r.branch(:, PF), [12.687; -30; -120; 240; 0; 120; 20; -62.3130; 82.3130; -42.687], 4, [t 'Pf']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 35.6504; 0; 0; 7.9756; 0; 0; 0; 0], 4, [t 'mu Pf']);
t_is(r.softlims.RATE_A.overload, zeros(10, 1), 12, [t 'softlims.overload']);
t_is(r.softlims.RATE_A.ovl_cost, zeros(10, 1), 12, [t 'softlims.ovl_cost']);
t_is(r.order.branch.status.on(r.order.int.softlims.RATE_A.idx), [3; 6; 8; 9; 10], 12, [t 'mu Pf']);
gen_order_check(mpc, r, t, mpopt1);

t = 'savecase(fname, mpc) : ';
fn = sprintf('softlims_savecase_test_%d', fix(1e8*rand));
savecase(fn, mpc);
mpc1 = loadcase(fn);
delete([fn '.m']);
t_ok(isfield(mpc1, 'softlims'), [t 'mpc.softlims']);
for lm = fieldnames(mpc.softlims).'
    lim = lm{:};
    t_ok(isfield(mpc1.softlims, lim), [t 'mpc.softlims.' lim])
    for f = fieldnames(mpc.softlims.(lim)).'
        field = f{:};
        if any(isinf(mpc.softlims.(lim).(field)))
            t_ok(isequal(mpc1.softlims.(lim).(field), mpc.softlims.(lim).(field)), [t 'mpc.softlims.' lim '.' field])
        else
            t_is(mpc1.softlims.(lim).(field), mpc.softlims.(lim).(field), 5, [t 'mpc.softlims.' lim '.' field])
        end
    end
end

t = 'savecase(fname, results)+gentype/fuel: ';
r.genfuel = {'hydro'; 'coal'; 'ng'; 'coal'};
r.gentype = {'HY'; 'ST'; 'GT'; 'ST'};
fn = sprintf('softlims_savecase_test_%d', fix(1e8*rand));
savecase(fn, r);
mpc1 = loadcase(fn);
delete([fn '.m']);
t_ok(isfield(mpc1, 'softlims'), [t 'results.softlims']);
for lm = fieldnames(r.softlims).'
    lim = lm{:};
    t_ok(isfield(mpc1.softlims, lim), [t 'results.softlims.' lim])
    for f = fieldnames(r.softlims.(lim)).'
        field = f{:};
        if any(isinf(r.softlims.(lim).(field)))
            t_ok(isequal(mpc1.softlims.(lim).(field), r.softlims.(lim).(field)), [t 'results.softlims.' lim '.' field])
        else
            t_is(mpc1.softlims.(lim).(field), r.softlims.(lim).(field), 5, [t 'results.softlims.' lim '.' field])
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

t = 'DC - softlims (violated) : ';
mpc = rmfield(mpc, 'softlims');
mpc.softlims.RATE_A.hl_mod = 'remove';
mpc.softlims.RATE_A.cost = 20;
for lim = {'VMIN', 'VMAX', 'ANGMIN', 'ANGMAX', 'PMIN', 'PMAX', 'QMIN', 'QMAX'}
    mpc.softlims.(lim{:}).hl_mod = 'none';
end
r = rundcopf(mpc, mpopt1);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9106.5059, 4, [t 'f']);
t_is(r.gen(:, PG), [241.3012         0   10.0000   63.6988].', 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -31.3012; -121.3012; 241.3012; 0; 120; 20; -63.6988; 83.6988; -41.3012], 4, [t 'Pf']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 20; 0; 0; 13.2992; 0; 0; 0; 0], 4, [t 'mu Pf']);
t_is(r.softlims.RATE_A.overload, [0; 0; 1.3011811; 0; 0; 0; 0; 0; 0; 0], 6, [t 'softlims.RATE_A.overload']);
t_is(r.softlims.RATE_A.ovl_cost, [0; 0; 26.023622; 0; 0; 0; 0; 0; 0; 0], 6, [t 'softlims.RATE_A.ovl_cost']);
gen_order_check(mpc, r, t, mpopt1);

t = 'savecase(fname, mpc) : ';
fn = sprintf('softlims_savecase_test_%d', fix(1e8*rand));
savecase(fn, mpc);
mpc1 = loadcase(fn);
delete([fn '.m']);
t_ok(isfield(mpc1, 'softlims'), [t 'mpc.softlims']);
for lm = fieldnames(mpc.softlims).'
    lim = lm{:};
    t_ok(isfield(mpc1.softlims, lim), [t 'mpc.softlims.' lim])
    for f = fieldnames(mpc.softlims.(lim)).'
        field = f{:};
        if any(isinf(mpc.softlims.(lim).(field)))
            t_ok(isequal(mpc1.softlims.(lim).(field), mpc.softlims.(lim).(field)), [t 'mpc.softlims.' lim '.' field])
        else
            t_is(mpc1.softlims.(lim).(field), mpc.softlims.(lim).(field), 5, [t 'mpc.softlims.' lim '.' field])
        end
    end
end

t = 'savecase(fname, results) : ';
fn = sprintf('softlims_savecase_test_%d', fix(1e8*rand));
savecase(fn, r);
mpc1 = loadcase(fn);
delete([fn '.m']);
t_ok(isfield(mpc1, 'softlims'), [t 'results.softlims']);
for lm = fieldnames(r.softlims).'
    lim = lm{:};
    t_ok(isfield(mpc1.softlims, lim), [t 'results.softlims.' lim])
    for f = fieldnames(r.softlims.(lim)).'
        field = f{:};
        if any(isinf(r.softlims.(lim).(field)))
            t_ok(isequal(mpc1.softlims.(lim).(field), r.softlims.(lim).(field)), [t 'results.softlims.' lim '.' field])
        else
            t_is(mpc1.softlims.(lim).(field), r.softlims.(lim).(field), 5, [t 'results.softlims.' lim '.' field])
        end
    end
end

%%-----  AC OPF (opf.flow_lim = 'S') -----
mpc = mpc0;
t = 'AC flow_lim=''S'' - hard limits : ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(~isfield(r.softlims.RATE_A, 'overload'), [t 'no softlims.RATE_A.overload']);
t_ok(~isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'no softlims.RATE_A.ovl_cost']);
t_is(r.f, 9521.385584, 4, [t 'f']);
t_is(r.gen(:, PG), [239.901164; 0 ; 17.312141; 66.456235], 4, [t 'Pg']);
t_is(r.branch(:, PF), [17.312141; -25.063413; -115.166831; 239.901164; 0; 119.964612; 18.530059; -66.456235; 84.941080; -42.201990], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-17.312141; 25.166831; 119.936552; -239.901164; 0; -118.530059; -18.484844; 66.456235; -82.798010; 42.375554], 4, [t 'Pt']);
t_is(r.branch(:, QF), [-2.721149; -2.842148; -16.241771; 21.275974; 0; -2.914100; -25.325409; -17.123405; 8.944864; -17.404939], 4, [t 'Qf']);
t_is(r.branch(:, QT), [2.884399; -13.758229; -3.901717; 6.815818; 0; -9.674591; 8.178541; 19.603112; -32.595061; -0.042250], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 29.242772; 0; 0; 8.561015; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
gen_order_check(mpc, r, t, mpopt);

t = 'AC flow_lim=''S'' - softlims(satisfied) : ';
mpc = toggle_softlims(mpc, 'on');
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9521.385584, 4, [t 'f']);
t_is(r.gen(:, PG), [239.901164; 0; 17.312141; 66.456235], 4, [t 'Pg']);
t_is(r.branch(:, PF), [17.312141; -25.063413; -115.166831; 239.901164; 0; 119.964612; 18.530059; -66.456235; 84.941080; -42.201990], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-17.312141; 25.166831; 119.936552; -239.901164; 0; -118.530059; -18.484844; 66.456235; -82.798010; 42.375554], 4, [t 'Pt']);
t_is(r.branch(:, QF), [-2.721149; -2.842148; -16.241771; 21.275974; 0; -2.914100; -25.325409; -17.123405; 8.944864; -17.404939], 4, [t 'Qf']);
t_is(r.branch(:, QT), [2.884399; -13.758229; -3.901717; 6.815818; 0; -9.674591; 8.178541; 19.603112; -32.595061; -0.042250], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 29.242772; 0; 0; 8.561015; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.RATE_A.overload, zeros(10, 1), 12, [t 'softlims.overload']);
t_is(r.softlims.RATE_A.ovl_cost, zeros(10, 1), 12, [t 'softlims.ovl_cost']);
t_is(r.order.branch.status.on(r.order.int.softlims.RATE_A.idx), [3; 6; 8; 9; 10], 12, [t 'softlims.RATE_A.idx']);
gen_order_check(mpc, r, t, mpopt);

t = 'AC flow_lim=''S'' - softlims (violated) : ';
mpc = rmfield(mpc, 'softlims');
mpc.softlims.RATE_A.hl_mod = 'remove';
mpc.softlims.RATE_A.cost = 20;
for lim = {'VMIN', 'VMAX', 'ANGMIN', 'ANGMAX', 'PMIN', 'PMAX', 'QMIN', 'QMAX'}
    mpc.softlims.(lim{:}).hl_mod = 'none';
end
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9486.265634, 4, [t 'f']);
t_is(r.gen(:, PG), [243.758931; 0; 10; 70.325987], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -28.621625; -118.761245; 243.758931; 0; 119.958848; 18.527937; -70.325987; 88.808549; -38.472172], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-10; 28.761245; 123.8; -243.758931; 0; -118.527937; -18.482561; 70.325987; -86.527828; 38.621625], 4, [t 'Pt']);
t_is(r.branch(:, QF), [3.364131; 0.552145; -12.847555; 19.471020; 0; -3.142423; -25.464124; -14.475416; 6.202740; -20.630111], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-3.306136; -17.152445; -6.346356; 9.488779; 0; -9.535876; 8.272676; 17.182536; -29.369889; 2.753991], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 20; 0; 0; 11.575125; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.RATE_A.overload, [0; 0; 3.962643; 0; 0; 0; 0; 0; 0; 0], 6, [t 'softlims.RATE_A.overload']);
t_is(r.softlims.RATE_A.ovl_cost, [0; 0; 79.252860; 0; 0; 0; 0; 0; 0; 0], 4, [t 'softlims.RATE_A.ovl_cost']);
gen_order_check(mpc, r, t, mpopt);

%% -----  AC OPF (opf.flow_lim = 'I') -----
mpopt = mpoption(mpopt, 'opf.flow_lim', 'I');
mpc = mpc0;
t = 'AC flow_lim=''I'' - hard limits : ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(~isfield(r.softlims.RATE_A, 'overload'), [t 'no softlims.RATE_A.overload']);
t_ok(~isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'no softlims.RATE_A.ovl_cost']);
t_is(r.f, 9178.300537, 4, [t 'f']);
t_is(r.gen(:, PG), [260.040338; 0; 10; 54.432302], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -32.676659; -122.910718; 260.040338; 0; 131.831065; 30.115930; -54.432302; 84.451894; -42.475047], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-10; 32.910718; 128.209273; -260.040338; 0; -130.115930; -30.019591; 54.432302; -82.524953; 42.676659], 4, [t 'Pt']);
t_is(r.branch(:, QF), [27.824349; 14.237742; 1.341441; 7.733221; 0; -4.887428; -29.462755; -4.743706; -7.797128; -31.778587], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-27.408203; -31.341441; -20.438656; 25.326084; 0; -5.537245; 12.540834; 6.294583; -18.221413; 13.170461], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 0; 0; 0; 20.114712; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
gen_order_check(mpc, r, t, mpopt);

t = 'AC flow_lim=''I'' - softlims(satisfied) : ';
mpc = toggle_softlims(mpc, 'on');
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9178.300537, 4, [t 'f']);
t_is(r.gen(:, PG), [260.040338; 0; 10; 54.432302], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -32.676659; -122.910718; 260.040338; 0; 131.831065; 30.115930; -54.432302; 84.451894; -42.475047], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-10; 32.910718; 128.209273; -260.040338; 0; -130.115930; -30.019591; 54.432302; -82.524953; 42.676659], 4, [t 'Pt']);
t_is(r.branch(:, QF), [27.824349; 14.237742; 1.341441; 7.733221; 0; -4.887428; -29.462755; -4.743706; -7.797128; -31.778587], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-27.408203; -31.341441; -20.438656; 25.326084; 0; -5.537245; 12.540834; 6.294583; -18.221413; 13.170461], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 0; 0; 0; 20.114712; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.RATE_A.overload, zeros(10, 1), 12, [t 'softlims.RATE_A.overload']);
t_is(r.softlims.RATE_A.ovl_cost, zeros(10, 1), 12, [t 'softlims.RATE_A.ovl_cost']);
t_is(r.order.branch.status.on(r.order.int.softlims.RATE_A.idx), [3; 6; 8; 9; 10], 12, [t 'softlims.RATE_A.idx']);
gen_order_check(mpc, r, t, mpopt);

t = 'AC flow_lim=''I'' - softlims (violated) : ';
mpc = rmfield(mpc, 'softlims');
mpc.softlims.RATE_A.hl_mod = 'remove';
mpc.softlims.RATE_A.cost = 15;
for lim = {'VMIN', 'VMAX', 'ANGMIN', 'ANGMAX', 'PMIN', 'PMAX', 'QMIN', 'QMAX'}
    mpc.softlims.(lim{:}).hl_mod = 'none';
end
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9146.413838, 4, [t 'f']);
t_is(r.gen(:, PG), [270; 0; 10; 44.862996], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -34.955827; -125.208395; 270; 0; 139.280487; 37.364683; -44.862996; 82.094229; -44.742338], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-10; 35.208395; 130.719513; -270; 0; -137.364683; -37.231233; 44.862996; -80.257662; 44.955827], 4, [t 'Pt']);
t_is(r.branch(:, QF), [25.262645; 13.437313; 0.303801; 13.552057; 0; -3.645909; -29.947558; -7.234346; -6.145528; -29.826268], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-24.907470; -30.303801; -18.341454; 21.987363; 0; -5.052442; 13.379873; 8.309624; -20.173732; 11.470157], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 6.260861; 0; 0; 15; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.RATE_A.overload, [0; 0; 0; 0; 0; 6.792934; 0; 0; 0; 0], 6, [t 'softlims.RATE_A.overload']);
t_is(r.softlims.RATE_A.ovl_cost, [0; 0; 0; 0; 0; 101.894009; 0; 0; 0; 0], 4, [t 'softlims.RATE_A.ovl_cost']);
gen_order_check(mpc, r, t, mpopt);

%%-----  AC OPF (opf.flow_lim = '2') -----
mpopt = mpoption(mpopt, 'opf.flow_lim', '2');
mpc = mpc0;
t = 'AC flow_lim=''2'' - hard limits : ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(~isfield(r.softlims.RATE_A, 'overload'), [t 'no softlims.RATE_A.overload']);
t_ok(~isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'no softlims.RATE_A.ovl_cost']);
t_is(r.f, 9513.613051, 4, [t 'f']);
t_is(r.gen(:, PG), [240; 0; 17.759246; 65.641269], 4, [t 'Pg']);
t_is(r.branch(:, PF), [17.759246; -25.216499; -115.346644; 240; 0; 120; 18.577311; -65.641269; 84.170657; -42.784248], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-17.759246; 25.346644; 120; -240; 0; -118.577311; -18.529387; 65.641269; -82.215752; 42.975745], 4, [t 'Pt']);
t_is(r.branch(:, QF), [15.646152; 6.577259; -6.087866; 8.021258; 0; -4.520662; -26.623563; -6.687069; -2.629290; -27.085546], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-15.370220; -23.912134; -15.548682; 20.069344; 0; -8.376437; 9.316359; 8.954090; -22.914454; 8.792961], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 29.248033; 0; 0; 8.417115; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
gen_order_check(mpc, r, t, mpopt);

t = 'AC flow_lim=''2'' - softlims(satisfied) : ';
mpc = toggle_softlims(mpc, 'on');
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9513.613051, 4, [t 'f']);
t_is(r.gen(:, PG), [240; 0; 17.759246; 65.641269], 4, [t 'Pg']);
t_is(r.branch(:, PF), [17.759246; -25.216499; -115.346644; 240; 0; 120; 18.577311; -65.641269; 84.170657; -42.784248], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-17.759246; 25.346644; 120; -240; 0; -118.577311; -18.529387; 65.641269; -82.215752; 42.975745], 4, [t 'Pt']);
t_is(r.branch(:, QF), [15.646152; 6.577259; -6.087866; 8.021258; 0; -4.520662; -26.623563; -6.687069; -2.629290; -27.085546], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-15.370220; -23.912134; -15.548682; 20.069344; 0; -8.376437; 9.316359; 8.954090; -22.914454; 8.792961], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 29.248033; 0; 0; 8.417115; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.RATE_A.overload, zeros(10, 1), 12, [t 'softlims.RATE_A.overload']);
t_is(r.softlims.RATE_A.ovl_cost, zeros(10, 1), 12, [t 'softlims.RATE_A.ovl_cost']);
t_is(r.order.branch.status.on(r.order.int.softlims.RATE_A.idx), [3; 6; 8; 9; 10], 12, [t 'softlims.RATE_A.idx']);
gen_order_check(mpc, r, t, mpopt);

t = 'AC flow_lim=''2'' - softlims (violated) : ';
mpc = rmfield(mpc, 'softlims');
mpc.softlims.RATE_A.hl_mod = 'remove';
mpc.softlims.RATE_A.cost = 20;
for lim = {'VMIN', 'VMAX', 'ANGMIN', 'ANGMAX', 'PMIN', 'PMAX', 'QMIN', 'QMAX'}
    mpc.softlims.(lim{:}).hl_mod = 'none';
end
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9476.854715, 4, [t 'f']);
t_is(r.gen(:, PG), [244.077740; 0; 10; 69.833910], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -28.933383; -119.111239; 244.077740; 0; 120; 18.578292; -69.833910; 88.362736; -38.762267], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-10; 29.111239; 124.077740; -244.077740; 0; -118.578292; -18.528826; 69.833910; -86.237733; 38.933383], 4, [t 'Pt']);
t_is(r.branch(:, QF), [22.204979; 10.325903; -2.434897; 5.980065; 0; -5.159694; -27.240912; -4.815545; -5.105567; -30.250169], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-21.917920; -27.565103; -17.970864; 23.130558; 0; -7.759088; 9.921112; 7.362539; -19.749831; 11.592017], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 20; 0; 0; 11.526783; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.RATE_A.overload, [0; 0; 4.07774; 0; 0; 0; 0; 0; 0; 0], 6, [t 'softlims.RATE_A.overload']);
t_is(r.softlims.RATE_A.ovl_cost, [0; 0; 81.554809; 0; 0; 0; 0; 0; 0; 0], 4, [t 'softlims.RATE_A.ovl_cost']);
gen_order_check(mpc, r, t, mpopt);

%%-----  AC OPF (opf.flow_lim = 'P') -----
mpopt = mpoption(mpopt, 'opf.flow_lim', 'P');
mpc = mpc0;
t = 'AC flow_lim=''P'' - hard limits : ';
t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(~isfield(r.softlims.RATE_A, 'overload'), [t 'no softlims.RATE_A.overload']);
t_ok(~isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'no softlims.RATE_A.ovl_cost']);
t_is(r.f, 9513.613051, 4, [t 'f']);
t_is(r.gen(:, PG), [240; 0; 17.759246; 65.641269], 4, [t 'Pg']);
t_is(r.branch(:, PF), [17.759246; -25.216499; -115.346644; 240; 0; 120; 18.577311; -65.641269; 84.170657; -42.784248], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-17.759246; 25.346644; 120; -240; 0; -118.577311; -18.529387; 65.641269; -82.215752; 42.975745], 4, [t 'Pt']);
t_is(r.branch(:, QF), [15.646152; 6.577259; -6.087866; 8.021258; 0; -4.520662; -26.623563; -6.687069; -2.629290; -27.085546], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-15.370220; -23.912134; -15.548682; 20.069344; 0; -8.376437; 9.316359; 8.954090; -22.914454; 8.792961], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 29.248033; 0; 0; 8.417115; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
gen_order_check(mpc, r, t, mpopt);

t = 'AC flow_lim=''P'' - softlims(satisfied) : ';
mpc = toggle_softlims(mpc, 'on');
t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9513.613051, 4, [t 'f']);
t_is(r.gen(:, PG), [240; 0; 17.759246; 65.641269], 4, [t 'Pg']);
t_is(r.branch(:, PF), [17.759246; -25.216499; -115.346644; 240; 0; 120; 18.577311; -65.641269; 84.170657; -42.784248], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-17.759246; 25.346644; 120; -240; 0; -118.577311; -18.529387; 65.641269; -82.215752; 42.975745], 4, [t 'Pt']);
t_is(r.branch(:, QF), [15.646152; 6.577259; -6.087866; 8.021258; 0; -4.520662; -26.623563; -6.687069; -2.629290; -27.085546], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-15.370220; -23.912134; -15.548682; 20.069344; 0; -8.376437; 9.316359; 8.954090; -22.914454; 8.792961], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 29.248033; 0; 0; 8.417115; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.RATE_A.overload, zeros(10, 1), 12, [t 'softlims.RATE_A.overload']);
t_is(r.softlims.RATE_A.ovl_cost, zeros(10, 1), 12, [t 'softlims.RATE_A.ovl_cost']);
t_is(r.order.branch.status.on(r.order.int.softlims.RATE_A.idx), [3; 6; 8; 9; 10], 12, [t 'softlims.RATE_A.idx']);
gen_order_check(mpc, r, t, mpopt);

t = 'AC flow_lim=''P'' - softlims (violated) : ';
mpc = rmfield(mpc, 'softlims');
mpc.softlims.RATE_A.hl_mod = 'remove';
mpc.softlims.RATE_A.cost = 20;
for lim = {'VMIN', 'VMAX', 'ANGMIN', 'ANGMAX', 'PMIN', 'PMAX', 'QMIN', 'QMAX'}
    mpc.softlims.(lim{:}).hl_mod = 'none';
end
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_ok(isfield(r.softlims.RATE_A, 'overload'), [t 'softlims.RATE_A.overload exists']);
t_ok(isfield(r.softlims.RATE_A, 'ovl_cost'), [t 'softlims.RATE_A.ovl_cost exists']);
t_is(r.f, 9476.854715, 4, [t 'f']);
t_is(r.gen(:, PG), [244.077740; 0; 10; 69.833910], 4, [t 'Pg']);
t_is(r.branch(:, PF), [10; -28.933383; -119.111239; 244.077740; 0; 120; 18.578292; -69.833910; 88.362736; -38.762267], 4, [t 'Pf']);
t_is(r.branch(:, PT), [-10; 29.111239; 124.077740; -244.077740; 0; -118.578292; -18.528826; 69.833910; -86.237733; 38.933383], 4, [t 'Pt']);
t_is(r.branch(:, QF), [22.204979; 10.325903; -2.434897; 5.980065; 0; -5.159694; -27.240912; -4.815545; -5.105567; -30.250169], 4, [t 'Qf']);
t_is(r.branch(:, QT), [-21.917920; -27.565103; -17.970864; 23.130558; 0; -7.759088; 9.921112; 7.362539; -19.749831; 11.592017], 4, [t 'Qt']);
t_is(r.branch(:, MU_SF)+r.branch(:, MU_ST), [0; 0; 20; 0; 0; 11.526783; 0; 0; 0; 0], 4, [t 'mu Sf+St']);
t_is(r.softlims.RATE_A.overload, [0; 0; 4.07774; 0; 0; 0; 0; 0; 0; 0], 6, [t 'softlims.RATE_A.overload']);
t_is(r.softlims.RATE_A.ovl_cost, [0; 0; 81.554809; 0; 0; 0; 0; 0; 0; 0], 4, [t 'softlims.RATE_A.ovl_cost']);
gen_order_check(mpc, r, t, mpopt);

if have_fcn('octave')
    warning(s1.state, file_in_path_warn_id);
end
warning(s2.state, sing_matrix_warn_id);

t_end;


function overload_loop(r, t, yn, lst, isDC)
%% r is the result structure
%% t is the test string
%% yn is whether to check if it is (1) or is not (0) a field
%% lst is a cell of limit names
%%     if yn == 1 i.e. checking if limits exist then this is a list of
%%     limits to check
%%     if yn == 0 i.e. checking that the limits do not exist, this is
%%     a list of limits to exclude

if nargin < 5
    isDC = 0;
    if nargin < 4
        lst = {};
    end
end

if isDC
    mpopt = mpoption('model', 'DC');
else
    mpopt = mpoption('model', 'AC');
end
lims_default = fieldnames(softlims_lim2mat(mpopt)).';

if ~yn
    lst = lims_default(~ismember(lims_default, lst));
else
    if isempty(lst)
        lst = lims_default;
    else
        lst = lims_default(ismember(lims_default, lst));
    end
end
for lm = lst
    lim = lm{1};
    if yn
        t_ok(isfield(r.softlims.(lim), 'overload'), [t sprintf('softlims.%s.overload exists', lim)]);
        t_ok(isfield(r.softlims.(lim), 'ovl_cost'), [t sprintf('softlims.%s.ovl_cost exists', lim)]);
    else
        t_ok(~isfield(r.softlims.(lim), 'overload'), [t sprintf('no softlims.%s.overload', lim)]);
        t_ok(~isfield(r.softlims.(lim), 'ovl_cost'), [t sprintf('no softlims.%s.ovl_cost', lim)]);
    end
end

function r = toggle_run_check(mpc, mpopt, t, on_off)
%% checks the softlims status of mpc, runs the opf, checks for success and
%% returns the results structure

if ~on_off
    t_ok(~toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 0');
else
    t_ok(toggle_softlims(mpc, 'status'), 'toggle_softlims(mpc, ''status'') == 1');
end
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);

function mu_cost_test(r, t, mpopt)
%% Tests whether the shadow price on violations matches the respective
%% softlim cost. The shadow price should equal the softlim cost when the
%% violation variable is used but the constraint is not binding. i.e. there is
%% a non-zero overload but the violation variable is not at its own limit.
%% the function therefore searches for any non-zero overload (mumask1)
%% that is also below the violation variables upper bound stored in s.ub
%% (mumask2). The cost of these variables is compared to the respective
%% shadow prices.

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% structures used to index into softlims in a loop
lims = softlims_lim2mat(mpopt);

vv = r.om.get_idx();
[x0, xmin, xmax] = r.om.params_var();

for lm = fieldnames(lims).'
    lim = lm{:};
    mat_name = lims.(lim).mat_name;
    if ~isfield(r.softlims, lim)
        continue
    end
    if strcmp(r.softlims.(lim).hl_mod, 'none')
        continue
    end
    s = r.softlims.(lim);

    varname = ['s_' lower(lim)];
    ub = xmax(vv.i1.(varname):vv.iN.(varname));
    switch lim
        case {'RATE_A', 'PMAX', 'PMIN', 'QMAX', 'QMIN'}
            ub = ub * r.baseMVA;
        case {'ANGMAX', 'ANGMIN'}
            ub = ub * 180/pi;
    end

    %% violation variable is non zero
    mumask1 = s.overload > 1e-6;
    %% ensure that overloads are not at the violation variable boundary since
    %% then the shadow price can take on any value again.
    mumask2 = true(size(mumask1));
    mumask2(s.idx) = s.overload(s.idx) < ub;
    
    mumask = mumask1 & mumask2;
    
    %% linear cost
    cst = zeros(size(mumask1));
    cst(s.idx) = s.cost;
    cst = cst(mumask); %keep only relevant entries
    
    if strcmp(lim, 'RATE_A')
        mu = sum(r.branch(mumask, MU_SF:MU_ST));
    else
        muidx = eval(['MU_' lim]);
        mu = r.(mat_name)(mumask, muidx);
    end
    
    t_is(mu,cst,4,[t 'mu ' lim '=cost'])
end

function gen_order_check(mpc, r, t, mpopt)
%% Make sure the ordering of generators in input case and result 
%% case is the same. Also make sure that any overloads are correctly
%% reflected in the difference between the generation value and
%% the respective limit

if nargin < 4
    isAC = 1;
else
    isAC = ~(~isempty(mpopt) && isfield(mpopt, 'model') && strcmp(mpopt.model, 'DC'));
end

[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

glist = {'P'};
if isAC
    t_is(mpc.gen(:,QMAX), r.gen(:,QMAX), 6, [t 'matching QMAX'])
    t_is(mpc.gen(:,QMIN), r.gen(:,QMIN), 6, [t 'matching QMIN'])
    glist = {glist{:}, 'Q'};
end
t_is(mpc.gen(:,PMAX), r.gen(:,PMAX), 6, [t 'matching PMAX'])
t_is(mpc.gen(:,PMIN), r.gen(:,PMIN), 6, [t 'matching PMIN'])
t_is(mpc.gen(:,GEN_STATUS), r.gen(:,GEN_STATUS), 6, [t 'matching gen status'])

on = r.gen(:,GEN_STATUS) > 0;
for gg = glist
    g = gg{:};
    if isfield(r, 'softlims')
        gidx = eval([g 'G']);
        gmax = eval([g 'MAX']);
        gmin = eval([g 'MIN']);
        if isfield(r.softlims.([g 'MAX']), 'overload')
            mask = (r.gen(:,gidx) - r.gen(:,gmax) > 0) & on;
            t_is(r.softlims.([g 'MAX']).overload(mask), r.gen(mask, gidx) - r.gen(mask,gmax), 4, [t g 'G-' g 'MAX= ' g 'MAX overload'])
        end
        if isfield(r.softlims.([g 'MIN']), 'overload')
            mask = (r.gen(:,gmin) - r.gen(:,gidx) > 0) & on;
            t_is(r.softlims.([g 'MIN']).overload(mask), r.gen(mask, gmin) - r.gen(mask,gidx), 4, [t g 'MIN-' g 'G= ' g 'MIN overload'])
        end
    end
end


%%-----  softlims_lim2mat  --------------------------------------------
function lim2mat = softlims_lim2mat(mpopt)

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

if nargin < 1
    mpopt = [];
end
if ~isempty(mpopt) && isfield(mpopt, 'model') && strcmp(mpopt.model, 'DC')
    lim2mat = struct(...
        'PMAX',   struct('mat_name', 'gen',    'col', PMAX,   'dir',  1 ), ...
        'PMIN',   struct('mat_name', 'gen',    'col', PMIN,   'dir', -1 ), ...
        'RATE_A', struct('mat_name', 'branch', 'col', RATE_A, 'dir',  1 ), ...
        'ANGMAX', struct('mat_name', 'branch', 'col', ANGMAX, 'dir',  1 ), ...
        'ANGMIN', struct('mat_name', 'branch', 'col', ANGMIN, 'dir', -1 ) ...
    );
else
    lim2mat = struct(...
        'VMAX',   struct('mat_name', 'bus',    'col', VMAX,   'dir',  1 ), ...
        'VMIN',   struct('mat_name', 'bus',    'col', VMIN,   'dir', -1 ), ...
        'PMAX',   struct('mat_name', 'gen',    'col', PMAX,   'dir',  1 ), ...
        'PMIN',   struct('mat_name', 'gen',    'col', PMIN,   'dir', -1 ), ...
        'QMAX',   struct('mat_name', 'gen',    'col', QMAX,   'dir',  1 ), ...
        'QMIN',   struct('mat_name', 'gen',    'col', QMIN,   'dir', -1 ), ...
        'RATE_A', struct('mat_name', 'branch', 'col', RATE_A, 'dir',  1 ), ...
        'ANGMAX', struct('mat_name', 'branch', 'col', ANGMAX, 'dir',  1 ), ...
        'ANGMIN', struct('mat_name', 'branch', 'col', ANGMIN, 'dir', -1 ) ...
    );
end
