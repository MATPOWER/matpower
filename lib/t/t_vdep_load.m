function t_vdep_load(quiet)
%T_VDEP_LOAD    Test voltage dependent ZIP load model for PF, CPF, OPF.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

pfalgs = {'NR', 'FDXB', 'FDBX'};
cpfparm = {'natural', 'arc-len', 'pseudo arc-len'};

t_begin(33+27*length(pfalgs)+39*2*length(cpfparm), quiet);

if quiet
    verbose = 0;
else
    verbose = 0;
end
verbose = 0;

casefile = 't_case9_opf';
mpopt = mpoption;
mpopt = mpoption(mpopt, 'out.gen', 1);
mpopt = mpoption(mpopt, 'pf.tol', 1e-9);
mpopt = mpoption(mpopt, 'opf.ac.solver', 'MIPS');
mpopt = mpoption(mpopt, 'opf.violation', 1e-6, 'mips.gradtol', 1e-9, ...
        'mips.comptol', 1e-8, 'mips.costtol', 1e-8);
mpopt = mpoption(mpopt, 'verbose', verbose);
if ~verbose
    mpopt = mpoption(mpopt, 'out.all', 0);
end

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% load base case file
mpc = loadcase(casefile);
nb = size(mpc.bus, 1);

t = 'PF - base case : ';
mpopt0 = mpopt;
r0 = runpf(mpc, mpopt);
t_ok(r0.success, [t 'success']);

for k = 1:length(pfalgs)
    mpopt = mpoption(mpopt0, 'pf.alg', pfalgs{k});

    t = sprintf('PF (%s) - constant power only : ', pfalgs{k});
    zip_w = [1 0 0];
    mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads', struct('pw', zip_w, 'qw', zip_w));
    r = runpf(mpc, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.bus, r0.bus, 6, [t 'bus']);
    t_is(r.gen, r0.gen, 6, [t 'gen']);
    t_is(r.branch, r0.branch, 6, [t 'branch']);

    t = sprintf('PF (%s) - constant current only : ', pfalgs{k});
    zip_w = [0 1 0];
    mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads', struct('pw', zip_w, 'qw', zip_w));
    r = runpf(mpc, mpopt);
    [tPd, tQd] = total_load(r, 'bus', [], mpopt);
    % [tPd, tQd] = total_load(r, 'bus');
    mpc1 = mpc;
    mpc1.bus(:, PD) = mpc1.bus(:, PD) .* r.bus(:, VM);
    mpc1.bus(:, QD) = mpc1.bus(:, QD) .* r.bus(:, VM);
    r1 = runpf(mpc1, mpopt0);
    t_ok(r1.success, [t 'success']);
    t_ok(r.success, [t 'success']);
    t_is(r.bus(:, [VM VA]), r1.bus(:, [VM VA]), 6, [t 'bus voltages']);
    t_is(tPd, r1.bus(:, PD), 6, [t 'bus P loads']);
    t_is(tPd, [0; 0; 0; 0; 87.937521; 0; 98.674071; 0; 120.041555], 6, [t 'bus P loads']);
    t_is(tQd, r1.bus(:, QD), 6, [t 'bus Q loads']);
    t_is(tQd, [0; 0; 0; 0; 29.312507; 0; 34.535925; 0; 48.016622], 6, [t 'bus Q loads']);
    t_is(r.gen(:, [PG QG]), r1.gen(:, [PG QG]), 6, [t 'gen dispatches']);
    t_is(r.branch, r1.branch, 6, [t 'branch']);

    t = sprintf('PF (%s) - constant impedance only : ', pfalgs{k});
    mpc1 = mpc;
    mpc1.bus(:, GS) =  mpc1.bus(:, PD);
    mpc1.bus(:, BS) = -mpc1.bus(:, QD);
    mpc1.bus(:, PD) = 0;
    mpc1.bus(:, QD) = 0;
    r1 = runpf(mpc1, mpopt0);
    t_ok(r1.success, [t 'success']);

    zip_w = [0 0 1];
    mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads', struct('pw', zip_w, 'qw', []));
    r = runpf(mpc, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.bus(:, [VM VA]), r1.bus(:, [VM VA]), 6, [t 'bus voltages']);
    t_is(r.gen(:, [PG QG]), r1.gen(:, [PG QG]), 6, [t 'gen dispatches']);
    t_is(r.branch, r1.branch, 6, [t 'branch']);

    t = sprintf('PF (%s) - combo ZIP loads : ', pfalgs{k});
    zip_w = [0.1 0.4 0.5];
    mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads.pw', zip_w);
    r = runpf(mpc, mpopt);
    [tPd, tQd] = total_load(r, 'bus', [], mpopt);
    mpc1 = mpc;
    Vm = r.bus(:, VM);
    scale = [Vm.^0 Vm Vm.^2] * zip_w';
    mpc1.bus(:, PD) = mpc1.bus(:, PD) .* scale;
    mpc1.bus(:, QD) = mpc1.bus(:, QD) .* scale;
    r1 = runpf(mpc1, mpopt0);
    t_ok(r1.success, [t 'success']);

    t_ok(r.success, [t 'success']);
    t_is(r.bus(:, [VM VA]), r1.bus(:, [VM VA]), 6, [t 'bus voltages']);
    t_is(tPd, r1.bus(:, PD), 6, [t 'bus P loads']);
    t_is(tPd, [0; 0; 0; 0; 87.205063; 0; 98.205247; 0; 118.314510], 6, [t 'bus P loads']);
    t_is(tQd, r1.bus(:, QD), 6, [t 'bus Q loads']);
    t_is(tQd, [0; 0; 0; 0; 29.068354; 0; 34.371836; 0; 47.325804], 6, [t 'bus Q loads']);
    t_is(r.gen(:, [PG QG]), r1.gen(:, [PG QG]), 6, [t 'gen dispatches']);
    t_is(r.branch, r1.branch, 6, [t 'branch']);
end

t = 'OPF - base case : ';
% mpopt0 = mpoption(mpopt0, 'verbose', 2);
mpopt = mpopt0;
r0 = runopf(mpc, mpopt);
t_ok(r0.success, [t 'success']);

t = 'OPF - constant power only : ';
zip_w = [1 0 0];
mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads.pw', zip_w);
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.f, r0.f, 6, [t 'f']);
t_is(r.bus, r0.bus, 6, [t 'bus']);
t_is(r.gen, r0.gen, 6, [t 'gen']);
t_is(r.branch, r0.branch, 6, [t 'branch']);

t = 'OPF - constant current only : ';
zip_w = [0 1 0];
mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads.pw', zip_w);
r = runopf(mpc, mpopt);
[tPd, tQd] = total_load(r, 'bus', [], mpopt);
% [tPd, tQd] = total_load(r, 'bus');
mpc1 = mpc;
mpc1.bus(:, PD) = mpc1.bus(:, PD) .* r.bus(:, VM);
mpc1.bus(:, QD) = mpc1.bus(:, QD) .* r.bus(:, VM);
k = find(mpc1.bus(:, PD));  %% buses with non-zero loads
mpc1.bus(k, VMAX) = r.bus(k, VM);
mpc1.bus(k, VMIN) = r.bus(k, VM);
r1 = runopf(mpc1, mpopt0);
t_ok(r1.success, [t 'success']);
t_ok(r.success, [t 'success']);
t_is(r.f, r1.f, 4, [t 'f']);
t_is(r.bus(:, [VM VA]), r1.bus(:, [VM VA]), 6, [t 'bus voltages']);
t_is(tPd, r1.bus(:, PD), 6, [t 'bus P loads']);
t_is(tPd, [0; 0; 0; 0; 81.708978; 0; 90.914185; 0; 112.500000], 6, [t 'bus P loads']);
t_is(tQd, r1.bus(:, QD), 6, [t 'bus Q loads']);
t_is(tQd, [0; 0; 0; 0; 27.236326; 0; 31.819965; 0; 45.000000], 6, [t 'bus Q loads']);
t_is(r.gen(:, [PG QG]), r1.gen(:, [PG QG]), 6, [t 'gen dispatches']);
t_is(r.branch(:, 1:ANGMAX), r1.branch(:, 1:ANGMAX), 6, [t 'branch']);

t = 'OPF - constant impedance only : ';
mpc1 = mpc;
mpc1.bus(:, GS) =  mpc1.bus(:, PD);
mpc1.bus(:, BS) = -mpc1.bus(:, QD);
mpc1.bus(:, PD) = 0;
mpc1.bus(:, QD) = 0;
r1 = runopf(mpc1, mpopt0);
t_ok(r1.success, [t 'success']);

zip_w = [0 0 1];
mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads.pw', zip_w);
r = runopf(mpc, mpopt);
t_ok(r.success, [t 'success']);
t_is(r.f, r1.f, 6, [t 'f']);
t_is(r.bus(:, [VM VA]), r1.bus(:, [VM VA]), 6, [t 'bus voltages']);
t_is(r.gen(:, [PG QG]), r1.gen(:, [PG QG]), 6, [t 'gen dispatches']);
t_is(r.branch, r1.branch, 6, [t 'branch']);

t = 'OPF - combo ZIP loads : ';
zip_w = [0.1 0.4 0.5];
mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads.pw', zip_w);
r = runopf(mpc, mpopt);
[tPd, tQd] = total_load(r, 'bus', [], mpopt);
mpc1 = mpc;
Vm = r.bus(:, VM);
scale = [Vm.^0 Vm Vm.^2] * zip_w';
mpc1.bus(:, PD) = mpc1.bus(:, PD) .* scale;
mpc1.bus(:, QD) = mpc1.bus(:, QD) .* scale;
k = find(mpc1.bus(:, PD));  %% buses with non-zero loads
mpc1.bus(k, VMAX) = r.bus(k, VM);
mpc1.bus(k, VMIN) = r.bus(k, VM);
r1 = runopf(mpc1, mpopt0);
t_ok(r1.success, [t 'success']);

t_ok(r.success, [t 'success']);
t_is(r.f, r1.f, 4, [t 'f']);
t_is(r.bus(:, [VM VA]), r1.bus(:, [VM VA]), 6, [t 'bus voltages']);
t_is(tPd, r1.bus(:, PD), 6, [t 'bus P loads']);
t_is(tPd, [0; 0; 0; 0; 78.909380; 0; 87.300375; 0; 108.125000], 6, [t 'bus P loads']);
t_is(tQd, r1.bus(:, QD), 6, [t 'bus Q loads']);
t_is(tQd, [0; 0; 0; 0; 26.303127; 0; 30.555131; 0; 43.250000], 6, [t 'bus Q loads']);
t_is(r.gen(:, [PG QG]), r1.gen(:, [PG QG]), 6, [t 'gen dispatches']);
t_is(r.branch(:, 1:ANGMAX), r1.branch(:, 1:ANGMAX), 6, [t 'branch']);


mpopt =  mpopt0;
    plot_nose_curve = 0;
    mpopt = mpoption(mpopt, 'cpf.stop_at', 1);
    mpopt = mpoption(mpopt, 'cpf.step', 0.05);
    %mpopt = mpoption(mpopt, 'cpf.adapt_step', 1);
    %mpopt = mpoption(mpopt, 'cpf.adapt_step_tol', 2e-5);
    mpopt = mpoption(mpopt, 'cpf.plot.level', plot_nose_curve);
    %mpopt = mpoption(mpopt, 'cpf.plot.bus', 9);
mpopt0 = mpopt;

    %% set up base and target cases
    mpcb = loadcase(casefile);
    mpct = mpcb;
    factor = 2.5 * 0.7;
    mpct.gen(:, [PG QG]) = mpct.gen(:, [PG QG]) * factor;
    mpct.bus(:, [PD QD]) = mpct.bus(:, [PD QD]) * factor;

    iterations = [
        20 20 20 20;
        21 22 22 22;
        21 22 22 22;
    ];

for k = 1:length(cpfparm)
    mpopt0 = mpoption(mpopt0, 'cpf.parameterization', k);

    t = sprintf('CPF (%s) - base case : ', cpfparm{k});
    r0 = runcpf(mpcb, mpct, mpopt0);
    t_ok(r0.success, [t 'success']);
    t_is(r0.cpf.iterations, iterations(k, 1), 12, [t 'iterations']);
    t_is(r0.cpf.max_lam, 1, 12, [t 'max_lam']);

    mpopt = mpopt0;
    t = sprintf('CPF (%s) - constant power only : ', cpfparm{k});
    zip_w = [1 0 0];
    mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads.pw', zip_w);
    r = runcpf(mpcb, mpct, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, r0.cpf.iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, r0.cpf.max_lam, 12, [t 'max_lam']);
    t_is(r.bus, r0.bus, 6, [t 'bus']);
    t_is(r.gen, r0.gen, 6, [t 'gen']);
    t_is(r.branch, r0.branch, 6, [t 'branch']);

    t = sprintf('CPF (%s) - constant current only : ', cpfparm{k});
    zip_w = [0 1 0];
    mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads.pw', zip_w);
    rb = runpf(mpcb, mpopt);
    r = runcpf(mpcb, mpct, mpopt);
    [tPd, tQd] = total_load(r, 'bus', [], mpopt);
    mpc1b = mpcb;
    mpc1b.bus(:, PD) = mpc1b.bus(:, PD) .* rb.bus(:, VM);
    mpc1b.bus(:, QD) = mpc1b.bus(:, QD) .* rb.bus(:, VM);
    mpc1t = mpct;
    mpc1t.bus(:, PD) = mpc1t.bus(:, PD) .* r.bus(:, VM);
    mpc1t.bus(:, QD) = mpc1t.bus(:, QD) .* r.bus(:, VM);
    r1 = runcpf(mpc1b, mpc1t, mpopt0);
    t_ok(r1.success, [t 'success']);
    t_is(r1.cpf.iterations, iterations(k,2), 12, [t 'iterations']);
    t_ok(r.success, [t 'success']);
    t_is(r.bus(:, [VM VA]), r1.bus(:, [VM VA]), 6, [t 'bus voltages']);
    t_is(tPd, r1.bus(:, PD), 6, [t 'bus P loads']);
    t_is(tPd, [0; 0; 0; 0; 143.494121; 0; 163.580859; 0; 191.911932], 6, [t 'bus P loads']);
    t_is(tQd, r1.bus(:, QD), 6, [t 'bus Q loads']);
    t_is(tQd, [0; 0; 0; 0; 47.831374; 0; 57.253301; 0; 76.764773], 6, [t 'bus Q loads']);
    t_is(r.gen(:, [PG QG]), r1.gen(:, [PG QG]), 6, [t 'gen dispatches']);
    t_is(r.branch, r1.branch, 6, [t 'branch']);

    t = sprintf('CPF (%s) - constant impedance only : ', cpfparm{k});
    zip_w = [0 0 1];
    mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads.pw', zip_w);
    rb = runpf(mpcb, mpopt);
    r = runcpf(mpcb, mpct, mpopt);
    [tPd, tQd] = total_load(r, 'bus', [], mpopt);
    mpc1b = mpcb;
    mpc1b.bus(:, PD) = mpc1b.bus(:, PD) .* rb.bus(:, VM).^2;
    mpc1b.bus(:, QD) = mpc1b.bus(:, QD) .* rb.bus(:, VM).^2;
    mpc1t = mpct;
    mpc1t.bus(:, PD) = mpc1t.bus(:, PD) .* r.bus(:, VM).^2;
    mpc1t.bus(:, QD) = mpc1t.bus(:, QD) .* r.bus(:, VM).^2;
    r1 = runcpf(mpc1b, mpc1t, mpopt0);
    t_ok(r1.success, [t 'success']);
    t_is(r1.cpf.iterations, iterations(k,3), 12, [t 'iterations']);
    t_ok(r.success, [t 'success']);
    t_is(r.bus(:, [VM VA]), r1.bus(:, [VM VA]), 6, [t 'bus voltages']);
    t_is(tPd, r1.bus(:, PD), 6, [t 'bus P loads']);
    t_is(tPd, [0; 0; 0; 0; 132.977098; 0; 154.969557; 0; 172.896177], 6, [t 'bus P loads']);
    t_is(tQd, r1.bus(:, QD), 6, [t 'bus Q loads']);
    t_is(tQd, [0; 0; 0; 0; 44.325699; 0; 54.239345; 0; 69.158471], 6, [t 'bus Q loads']);
    t_is(r.gen(:, [PG QG]), r1.gen(:, [PG QG]), 6, [t 'gen dispatches']);
    t_is(r.branch, r1.branch, 6, [t 'branch']);

    t = sprintf('CPF (%s) - combo ZIP loads : ', cpfparm{k});
    zip_w = [0.1 0.4 0.5];
    mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads.pw', zip_w);
    rb = runpf(mpcb, mpopt);
    r = runcpf(mpcb, mpct, mpopt);
    [tPd, tQd] = total_load(r, 'bus', [], mpopt);
    Vm = abs(rb.bus(:, VM));
    scaleb = [Vm.^0 Vm Vm.^2] * zip_w';
    Vm = abs(r.bus(:, VM));
    scalet = [Vm.^0 Vm Vm.^2] * zip_w';
    mpc1b = mpcb;
    mpc1b.bus(:, PD) = mpc1b.bus(:, PD) .* scaleb;
    mpc1b.bus(:, QD) = mpc1b.bus(:, QD) .* scaleb;
    mpc1t = mpct;
    mpc1t.bus(:, PD) = mpc1t.bus(:, PD) .* scalet;
    mpc1t.bus(:, QD) = mpc1t.bus(:, QD) .* scalet;
    r1 = runcpf(mpc1b, mpc1t, mpopt0);
    t_ok(r1.success, [t 'success']);
    t_is(r1.cpf.iterations, iterations(k,4), 12, [t 'iterations']);
    t_ok(r.success, [t 'success']);
    t_is(r.bus(:, [VM VA]), r1.bus(:, [VM VA]), 6, [t 'bus voltages']);
    t_is(tPd, r1.bus(:, PD), 6, [t 'bus P loads']);
    t_is(tPd, [0; 0; 0; 0; 139.202168; 0; 160.005294; 0; 184.209518], 6, [t 'bus P loads']);
    t_is(tQd, r1.bus(:, QD), 6, [t 'bus Q loads']);
    t_is(tQd, [0; 0; 0; 0; 46.400723; 0; 56.001853; 0; 73.683807], 6, [t 'bus Q loads']);
    t_is(r.gen(:, [PG QG]), r1.gen(:, [PG QG]), 6, [t 'gen dispatches']);
    t_is(r.branch, r1.branch, 6, [t 'branch']);
end


%% set up base and target cases
    mpcb = loadcase(casefile);
    mpct = mpcb;
    factor = 2.5 * 0.6;
    mpct.gen(:, [PG QG]) = mpct.gen(:, [PG QG]) * factor;
    mpct.bus(:, [PD QD]) = mpct.bus(:, [PD QD]) * factor;
    mpct.gen(:, QMAX) = 75;
    mpcb.gen(:, QMAX) = 75;

    %% with Q limits enforced
    iterations = [
        22 20 20 20;
        22 21 21 21;
        22 21 21 21;
    ];

for k = 1:length(cpfparm)
    mpopt0 = mpoption(mpopt0, 'cpf.parameterization', k, 'cpf.enforce_q_lims', 1);

    t = sprintf('CPF (%s) w/Q lims - base case : ', cpfparm{k});
    r0 = runcpf(mpcb, mpct, mpopt0);
    t_ok(r0.success, [t 'success']);
    t_is(r0.cpf.iterations, iterations(k, 1), 12, [t 'iterations']);
    t_is(r0.cpf.max_lam, 1, 12, [t 'max_lam']);

    mpopt = mpopt0;
    t = sprintf('CPF (%s) w/Q lims - constant power only : ', cpfparm{k});
    zip_w = [1 0 0];
    mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads.pw', zip_w);
    r = runcpf(mpcb, mpct, mpopt);
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, r0.cpf.iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, r0.cpf.max_lam, 12, [t 'max_lam']);
    t_is(r.bus, r0.bus, 6, [t 'bus']);
    t_is(r.gen, r0.gen, 6, [t 'gen']);
    t_is(r.branch, r0.branch, 6, [t 'branch']);

    t = sprintf('CPF (%s) w/Q lims - constant current only : ', cpfparm{k});
    zip_w = [0 1 0];
    mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads.pw', zip_w);
    rb = runpf(mpcb, mpopt);
    r = runcpf(mpcb, mpct, mpopt);
    [tPd, tQd] = total_load(r, 'bus', [], mpopt);
    mpc1b = mpcb;
    mpc1b.bus(:, PD) = mpc1b.bus(:, PD) .* rb.bus(:, VM);
    mpc1b.bus(:, QD) = mpc1b.bus(:, QD) .* rb.bus(:, VM);
    mpc1t = mpct;
    mpc1t.bus(:, PD) = mpc1t.bus(:, PD) .* r.bus(:, VM);
    mpc1t.bus(:, QD) = mpc1t.bus(:, QD) .* r.bus(:, VM);
    r1 = runcpf(mpc1b, mpc1t, mpopt0);
    t_ok(r1.success, [t 'success']);
    t_is(r1.cpf.iterations, iterations(k,2), 12, [t 'iterations']);
    t_ok(r.success, [t 'success']);
    t_is(r.bus(:, [VM VA]), r1.bus(:, [VM VA]), 6, [t 'bus voltages']);
    t_is(tPd, r1.bus(:, PD), 6, [t 'bus P loads']);
    t_is(tPd, [0; 0; 0; 0; 126.246737; 0; 143.114682; 0; 170.231775], 6, [t 'bus P loads']);
    t_is(tQd, r1.bus(:, QD), 6, [t 'bus Q loads']);
    t_is(tQd, [0; 0; 0; 0; 42.082245; 0; 50.090138; 0; 68.092710], 6, [t 'bus Q loads']);
    t_is(r.gen(:, [PG QG]), r1.gen(:, [PG QG]), 6, [t 'gen dispatches']);
    t_is(r.branch, r1.branch, 6, [t 'branch']);

    t = sprintf('CPF (%s) w/Q lims - constant impedance only : ', cpfparm{k});
    zip_w = [0 0 1];
    mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads.pw', zip_w);
    rb = runpf(mpcb, mpopt);
    r = runcpf(mpcb, mpct, mpopt);
    [tPd, tQd] = total_load(r, 'bus', [], mpopt);
    mpc1b = mpcb;
    mpc1b.bus(:, PD) = mpc1b.bus(:, PD) .* rb.bus(:, VM).^2;
    mpc1b.bus(:, QD) = mpc1b.bus(:, QD) .* rb.bus(:, VM).^2;
    mpc1t = mpct;
    mpc1t.bus(:, PD) = mpc1t.bus(:, PD) .* r.bus(:, VM).^2;
    mpc1t.bus(:, QD) = mpc1t.bus(:, QD) .* r.bus(:, VM).^2;
    r1 = runcpf(mpc1b, mpc1t, mpopt0);
    t_ok(r1.success, [t 'success']);
    t_is(r1.cpf.iterations, iterations(k,3), 12, [t 'iterations']);
    t_ok(r.success, [t 'success']);
    t_is(r.bus(:, [VM VA]), r1.bus(:, [VM VA]), 6, [t 'bus voltages']);
    t_is(tPd, r1.bus(:, PD), 6, [t 'bus P loads']);
    t_is(tPd, [0; 0; 0; 0; 119.354069; 0; 137.673893; 0; 157.193899], 6, [t 'bus P loads']);
    t_is(tQd, r1.bus(:, QD), 6, [t 'bus Q loads']);
    t_is(tQd, [0; 0; 0; 0; 39.784689; 0; 48.185862; 0; 62.877559], 6, [t 'bus Q loads']);
    t_is(r.gen(:, [PG QG]), r1.gen(:, [PG QG]), 6, [t 'gen dispatches']);
    t_is(r.branch, r1.branch, 6, [t 'branch']);

    t = sprintf('CPF (%s) w/Q lims - combo ZIP loads : ', cpfparm{k});
    zip_w = [0.1 0.4 0.5];
    mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads.pw', zip_w);
    rb = runpf(mpcb, mpopt);
    r = runcpf(mpcb, mpct, mpopt);
    [tPd, tQd] = total_load(r, 'bus', [], mpopt);
    Vm = abs(rb.bus(:, VM));
    scaleb = [Vm.^0 Vm Vm.^2] * zip_w';
    Vm = abs(r.bus(:, VM));
    scalet = [Vm.^0 Vm Vm.^2] * zip_w';
    mpc1b = mpcb;
    mpc1b.bus(:, PD) = mpc1b.bus(:, PD) .* scaleb;
    mpc1b.bus(:, QD) = mpc1b.bus(:, QD) .* scaleb;
    mpc1t = mpct;
    mpc1t.bus(:, PD) = mpc1t.bus(:, PD) .* scalet;
    mpc1t.bus(:, QD) = mpc1t.bus(:, QD) .* scalet;
    r1 = runcpf(mpc1b, mpc1t, mpopt0);
    t_ok(r1.success, [t 'success']);
    t_is(r1.cpf.iterations, iterations(k,4), 12, [t 'iterations']);
    t_ok(r.success, [t 'success']);
    t_is(r.bus(:, [VM VA]), r1.bus(:, [VM VA]), 6, [t 'bus voltages']);
    t_is(tPd, r1.bus(:, PD), 6, [t 'bus P loads']);
    t_is(tPd, [0; 0; 0; 0; 123.418124; 0; 140.853612; 0; 164.913503], 6, [t 'bus P loads']);
    t_is(tQd, r1.bus(:, QD), 6, [t 'bus Q loads']);
    t_is(tQd, [0; 0; 0; 0; 41.1393746; 0; 49.298764; 0; 65.9654012], 6, [t 'bus Q loads']);
    t_is(r.gen(:, [PG QG]), r1.gen(:, [PG QG]), 6, [t 'gen dispatches']);
    t_is(r.branch, r1.branch, 6, [t 'branch']);
end

t_end;
