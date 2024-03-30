function res = t_cpf(quiet)
% t_cpf - Tests for legacy continuation power flow.

%   MATPOWER
%   Copyright (c) 2013-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

num_tests = 415;
t_begin(num_tests, quiet);

if have_feature('matlab', 'vnum') < 7.001
    t_skip(num_tests, 'RUNCPF requires cellfun() construct not available before MATLAB 7.1');
else
    plot_nose_curve = 0;
    verbose = 0;

    casefile = 't_case9_pfv2';
    if have_feature('octave')
        if have_feature('octave', 'vnum') >= 4
            file_in_path_warn_id = 'Octave:data-file-in-path';
        else
            file_in_path_warn_id = 'Octave:load-file-in-path';
        end
        s1 = warning('query', file_in_path_warn_id);
        warning('off', file_in_path_warn_id);
    end
    mpopt = mpoption('out.all', 0, 'verbose', verbose);
    %mpopt = mpoption(mpopt, 'cpf.stop_at', 'FULL', );
    mpopt = mpoption(mpopt, 'cpf.step', 0.02);
    %mpopt = mpoption(mpopt, 'cpf.adapt_step', 1);
    %mpopt = mpoption(mpopt, 'cpf.adapt_step_damping', 1);
    %mpopt = mpoption(mpopt, 'cpf.adapt_step_tol', 2e-5);
    mpopt = mpoption(mpopt, 'cpf.plot.level', plot_nose_curve);
    %mpopt = mpoption(mpopt, 'cpf.plot.bus', 9);
    %mpopt = mpoption(mpopt, 'pf.tol', 1e-10);
    %mpopt = mpoption(mpopt, 'verbose', 3);

    %% define named indices into bus, gen, branch matrices
    [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
        VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
    [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
        TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
        ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
    [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
        MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
        QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

    %% set up base and target cases
    mpcb = loadcase(casefile);
    %% add isolated bus to make sure int2ext works for V_hat, V
    mpcb.bus = [mpcb.bus(1:3, :); mpcb.bus(3, :); mpcb.bus(4:end, :)];
    mpcb.bus(4, BUS_I) = 50;
    mpcb.bus(4, BUS_TYPE) = NONE;
    mpcb.gen(1, QMAX) = 200;    %% decrease a Q lim
    % r = runpf(mpcb, mpopt);
    % mpcb.gen(1, [PG QG]) = r.gen(1, [PG QG]); %% solved values for slack gen
    mpct = mpcb;
    factor = 2.5;
    mpct.gen(:, [PG QG]) = mpct.gen(:, [PG QG]) * factor;
    mpct.bus(:, [PD QD]) = mpct.bus(:, [PD QD]) * factor;
    Pdb = mpcb.bus(:, PD);
    Pdt = mpct.bus(:, PD);
    Pd_lam = @(lam)(Pdb + lam * (Pdt - Pdb));
    pvg = [2;3];
    Pgb = mpcb.gen(pvg, PG);
    Pgt = mpct.gen(pvg, PG);
    Pg_lam = @(lam)(Pgb + lam * (Pgt - Pgb));

    %% run CPF
    t = 'base == target : ';
    r = runcpf(mpcb, mpcb, mpopt);
    iterations = 1;
    t_ok(r.success, [t 'success']);
    t_is(r.gen(pvg, PG), Pg_lam(0), 12, [t, 'Pg']);
    t_is(r.bus(:, PD), Pd_lam(0), 12, [t, 'Pd']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0, 12, [t 'max_lam']);
    t_is(size(r.cpf.V_hat), [10 1], 12, [t 'size(V_hat)']);
    t_is(size(r.cpf.V), [10 1], 12, [t 'size(V)']);
    t_is(size(r.cpf.lam_hat), [1 1], 12, [t 'size(lam_hat)']);
    t_is(size(r.cpf.lam), [1 1], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'Base case and target case have identical load and generation'), [t 'done_msg']);
    t_is(length(r.cpf.events), 0, 12, [t 'length(events) == 0']);

    t = 'CPF to lambda = 0.7 (natural) : ';
    mpopt = mpoption(mpopt, 'cpf.stop_at', 0.7, 'cpf.parameterization', 1);
    r = runcpf(mpcb, mpct, mpopt);
    iterations = 35;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0.7, 12, [t 'max_lam']);
    t_is(r.gen(pvg, PG), Pg_lam(0.7), 12, [t, 'Pg']);
    t_is(r.bus(:, PD), Pd_lam(0.7), 12, [t, 'Pd']);
    t_is(size(r.cpf.V_hat), [10 iterations+1], 12, [t 'size(V_hat)']);
    t_is(size(r.cpf.V), [10 iterations+1], 12, [t 'size(V)']);
    t_is(size(r.cpf.lam_hat), [1 iterations+1], 12, [t 'size(lam_hat)']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'Reached desired lambda 0.7'), [t 'done_msg']);
    t_is(length(r.cpf.events), 1, 12, [t 'length(events) == 1']);
    t_is(r.cpf.events(1).k, iterations, 12, [t 'events(1).k']);
    t_is(r.cpf.events(1).idx, 1, 12, [t 'events(1).idx']);
    t_str_match(r.cpf.events(1).name, 'TARGET_LAM', [t 'events(1).name']);

    t = 'CPF to lambda = 0.7 (arc length) : ';
    mpopt = mpoption(mpopt, 'cpf.stop_at', 0.7, 'cpf.parameterization', 2);
    r = runcpf(mpcb, mpct, mpopt);
    iterations = 41;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0.7, 12, [t 'max_lam']);
    t_is(r.gen(pvg, PG), Pg_lam(0.7), 12, [t, 'Pg']);
    t_is(r.bus(:, PD), Pd_lam(0.7), 12, [t, 'Pd']);
    t_is(size(r.cpf.V_hat), [10 iterations+1], 12, [t 'size(V_hat)']);
    t_is(size(r.cpf.V), [10 iterations+1], 12, [t 'size(V)']);
    t_is(size(r.cpf.lam_hat), [1 iterations+1], 12, [t 'size(lam_hat)']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'Reached desired lambda 0.7'), [t 'done_msg']);
    t_is(length(r.cpf.events), 1, 12, [t 'length(events) == 1']);
    t_is(r.cpf.events(1).k, iterations, 12, [t 'events(1).k']);
    t_is(r.cpf.events(1).idx, 1, 12, [t 'events(1).idx']);
    t_str_match(r.cpf.events(1).name, 'TARGET_LAM', [t 'events(1).name']);

    t = 'CPF to lambda = 0.7 (pseudo arc length) : ';
    mpopt = mpoption(mpopt, 'cpf.stop_at', 0.7, 'cpf.parameterization', 3);
    r = runcpf(mpcb, mpct, mpopt);
    iterations = 41;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0.7, 12, [t 'max_lam']);
    t_is(r.gen(pvg, PG), Pg_lam(0.7), 12, [t, 'Pg']);
    t_is(r.bus(:, PD), Pd_lam(0.7), 12, [t, 'Pd']);
    t_is(size(r.cpf.V_hat), [10 iterations+1], 12, [t 'size(V_hat)']);
    t_is(size(r.cpf.V), [10 iterations+1], 12, [t 'size(V)']);
    t_is(size(r.cpf.lam_hat), [1 iterations+1], 12, [t 'size(lam_hat)']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'Reached desired lambda 0.7'), [t 'done_msg']);
    t_is(length(r.cpf.events), 1, 12, [t 'length(events) == 1']);
    t_is(r.cpf.events(1).k, iterations, 12, [t 'events(1).k']);
    t_is(r.cpf.events(1).idx, 1, 12, [t 'events(1).idx']);
    t_str_match(r.cpf.events(1).name, 'TARGET_LAM', [t 'events(1).name']);

    t = 'CPF to nose pt (arc length) : ';
    mpopt = mpoption(mpopt, 'cpf.stop_at', 'NOSE', 'cpf.parameterization', 2);
    mpopt = mpoption(mpopt, 'cpf.adapt_step', 1);
    r = runcpf(mpcb, mpct, mpopt);
    iterations = 23;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0.99025, 3, [t 'max_lam']);
    t_is(r.gen(pvg, PG), Pg_lam(r.cpf.max_lam), 12, [t, 'Pg']);
    t_is(r.bus(:, PD), Pd_lam(r.cpf.max_lam), 12, [t, 'Pd']);
    t_is(size(r.cpf.V_hat), [10 iterations+1], 12, [t 'size(V_hat)']);
    t_is(size(r.cpf.V), [10 iterations+1], 12, [t 'size(V)']);
    t_is(size(r.cpf.lam_hat), [1 iterations+1], 12, [t 'size(lam_hat)']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'Reached steady state loading limit'), [t 'done_msg']);
    t_is(length(r.cpf.events), 1, 12, [t 'length(events) == 1']);
    t_is(r.cpf.events(1).k, iterations, 12, [t 'events(1).k']);
    t_is(r.cpf.events(1).idx, 1, 12, [t 'events(1).idx']);
    t_str_match(r.cpf.events(1).name, 'NOSE', [t 'events(1).name']);

    t = 'CPF to nose pt (pseudo arc length) : ';
    mpopt = mpoption(mpopt, 'cpf.stop_at', 'NOSE', 'cpf.parameterization', 3);
    mpopt = mpoption(mpopt, 'cpf.adapt_step', 1);
    r = runcpf(mpcb, mpct, mpopt);
    iterations = 23;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0.99025, 3, [t 'max_lam']);
    t_is(r.gen(pvg, PG), Pg_lam(r.cpf.max_lam), 12, [t, 'Pg']);
    t_is(r.bus(:, PD), Pd_lam(r.cpf.max_lam), 12, [t, 'Pd']);
    t_is(size(r.cpf.V_hat), [10 iterations+1], 12, [t 'size(V_hat)']);
    t_is(size(r.cpf.V), [10 iterations+1], 12, [t 'size(V)']);
    t_is(size(r.cpf.lam_hat), [1 iterations+1], 12, [t 'size(lam_hat)']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'Reached steady state loading limit'), [t 'done_msg']);
    t_is(length(r.cpf.events), 1, 12, [t 'length(events) == 1']);
    t_is(r.cpf.events(1).k, iterations, 12, [t 'events(1).k']);
    t_is(r.cpf.events(1).idx, 1, 12, [t 'events(1).idx']);
    t_str_match(r.cpf.events(1).name, 'NOSE', [t 'events(1).name']);

    t = 'CPF to nose pt (pseudo arc length) w/Q lims: ';
    mpopt_qlim = mpoption(mpopt, 'cpf.stop_at', 'NOSE', 'cpf.parameterization', 3,'cpf.enforce_q_lims',1);
    mpopt_qlim = mpoption(mpopt_qlim, 'cpf.adapt_step', 1);
    r = runcpf(mpcb, mpct, mpopt_qlim);
    iterations = 19;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0.795809, 6, [t 'max_lam']);
    ePg = Pg_lam(r.cpf.max_lam);
    t_is(r.gen(pvg(1), PG), ePg(1), 12, [t, 'Pg']);
    t_is(r.bus(:, PD), Pd_lam(r.cpf.max_lam), 12, [t, 'Pd']);
    t_is(size(r.cpf.V_hat), [10 iterations+1], 12, [t 'size(V_hat)']);
    t_is(size(r.cpf.V), [10 iterations+1], 12, [t 'size(V)']);
    t_is(size(r.cpf.lam_hat), [1 iterations+1], 12, [t 'size(lam_hat)']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'Reached steady state loading limit'), [t 'done_msg']);
    t_is(r.bus(:, BUS_TYPE), [1;3;2;4;1;1;1;1;1;1], 12, [t 'bus type']);
    ek = [12 iterations];
    eidx = [1 1];
    ename = {'QLIM', 'NOSE'};
    ne = length(ek);
    t_is(length(r.cpf.events), ne, 12, sprintf('%ssize(events) == %d', t, ne));
    for j = 1:ne
        t_is(r.cpf.events(j).k, ek(j), 12, sprintf('%sevents(%d).k == %d', t, j, ek(j)));
        t_is(r.cpf.events(j).idx, eidx(j), 12, sprintf('%sevents(%d).idx == %d', t, j, eidx(j)));
        t_str_match(r.cpf.events(j).name, ename{j}, sprintf('%sevents(%d).name = ''%s''', t, j, ename{j}));
    end
    
    t = 'CPF to nose pt (pseudo arc length) w/P lims: ';
    mpopt_plim = mpoption(mpopt, 'cpf.stop_at', 'NOSE', 'cpf.parameterization', 3,'cpf.enforce_p_lims',1);
    mpopt_plim = mpoption(mpopt_plim, 'cpf.adapt_step', 1);
    r = runcpf(mpcb, mpct, mpopt_plim);
    iterations = 19;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0.97975, 4, [t 'max_lam']);
    t_is(r.gen(pvg(2), PG), 300, 12, [t, 'Pg']);
    t_is(r.bus(:, PD), Pd_lam(r.cpf.max_lam), 12, [t, 'Pd']);
    t_is(size(r.cpf.V_hat), [10 iterations+1], 12, [t 'size(V_hat)']);
    t_is(size(r.cpf.V), [10 iterations+1], 12, [t 'size(V)']);
    t_is(size(r.cpf.lam_hat), [1 iterations+1], 12, [t 'size(lam_hat)']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_str_match(r.cpf.done_msg, 'All generators at Pmax', [t 'done_msg']);
    t_is(r.bus(:, BUS_TYPE), [2;2;3;4;1;1;1;1;1;1], 12, [t 'bus type']);
    ek = [7 13 iterations];
    eidx = [3 1 2];
    ename = {'PLIM', 'PLIM', 'PLIM'};
    ne = length(ek);
    t_is(length(r.cpf.events), ne, 12, sprintf('%ssize(events) == %d', t, ne));
    for j = 1:ne
        t_is(r.cpf.events(j).k, ek(j), 12, sprintf('%sevents(%d).k == %d', t, j, ek(j)));
        t_is(r.cpf.events(j).idx, eidx(j), 12, sprintf('%sevents(%d).idx == %d', t, j, eidx(j)));
        t_str_match(r.cpf.events(j).name, ename{j}, sprintf('%sevents(%d).name = ''%s''', t, j, ename{j}));
    end
    
    t = 'CPF to nose pt (pseudo arc length) w/flow lims: ';
    mpopt_flowlim = mpoption(mpopt, 'cpf.stop_at', 'NOSE', 'cpf.parameterization', 3,'cpf.enforce_flow_lims',1);
    mpopt_flowlim = mpoption(mpopt_flowlim, 'cpf.adapt_step', 1);
    r = runcpf(mpcb, mpct, mpopt_flowlim);
    iterations = 3;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0.110684, 6, [t 'max_lam']);
    t_is(r.gen(pvg, PG), Pg_lam(r.cpf.max_lam), 12, [t, 'Pg']);
    t_is(r.bus(:, PD), Pd_lam(r.cpf.max_lam), 12, [t, 'Pd']);
    t_is(size(r.cpf.V_hat), [10 iterations+1], 12, [t 'size(V_hat)']);
    t_is(size(r.cpf.V), [10 iterations+1], 12, [t 'size(V)']);
    t_is(size(r.cpf.lam_hat), [1 iterations+1], 12, [t 'size(lam_hat)']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'branch flow limit reached'), [t 'done_msg']);
    t_is(length(r.cpf.events), 1, 12, [t 'length(events) == 1']);
    t_is(r.cpf.events(1).k, iterations, 12, [t 'events(1).k']);
    t_is(r.cpf.events(1).idx, 5, 12, [t 'events(1).idx']);
    t_str_match(r.cpf.events(1).name, 'FLIM', [t 'events(1).name']);
    
    t = 'CPF to nose pt (pseudo arc length) w/violated flow lims: ';
    mpopt_flowlim = mpoption(mpopt, 'cpf.stop_at', 'NOSE', 'cpf.parameterization', 3,'cpf.enforce_flow_lims',1);
    mpopt_flowlim = mpoption(mpopt_flowlim, 'cpf.adapt_step', 1);
    mpcb1 = mpcb;
    mpcb1.branch(1, RATE_A) = 75;
    r = runcpf(mpcb1, mpct, mpopt_flowlim);
    iterations = 1;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0, 6, [t 'max_lam']);
    t_is(r.gen(pvg, PG), Pg_lam(r.cpf.max_lam), 12, [t, 'Pg']);
    t_is(r.bus(:, PD), Pd_lam(r.cpf.max_lam), 12, [t, 'Pd']);
    t_is(size(r.cpf.V_hat), [10 iterations], 12, [t 'size(V_hat)']);
    t_is(size(r.cpf.V), [10 iterations], 12, [t 'size(V)']);
    t_is(size(r.cpf.lam_hat), [1 iterations], 12, [t 'size(lam_hat)']);
    t_is(size(r.cpf.lam), [1 iterations], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'branch flow limit violated in base case'), [t 'done_msg']);
    t_is(length(r.cpf.events), 0, 12, [t 'length(events) == 0']);
    
    t = 'CPF to nose pt (pseudo arc length) w/V lims: ';
    mpopt_vlim = mpoption(mpopt, 'cpf.stop_at', 'NOSE', 'cpf.parameterization', 3,'cpf.enforce_v_lims',1);
    mpopt_vlim = mpoption(mpopt_vlim, 'cpf.adapt_step', 1);
    r = runcpf(mpcb, mpct, mpopt_vlim);
    iterations = 6;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0.316932, 6, [t 'max_lam']);
    t_is(r.gen(pvg, PG), Pg_lam(r.cpf.max_lam), 12, [t, 'Pg']);
    t_is(r.bus(:, PD), Pd_lam(r.cpf.max_lam), 12, [t, 'Pd']);
    t_is(size(r.cpf.V_hat), [10 iterations+1], 12, [t 'size(V_hat)']);
    t_is(size(r.cpf.V), [10 iterations+1], 12, [t 'size(V)']);
    t_is(size(r.cpf.lam_hat), [1 iterations+1], 12, [t 'size(lam_hat)']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'voltage magnitude limit reached'), [t 'done_msg']);
    t_is(length(r.cpf.events), 1, 12, [t 'length(events) == 1']);
    t_is(r.cpf.events(1).k, iterations, 12, [t 'events(1).k']);
    t_is(r.cpf.events(1).idx, 9, 12, [t 'events(1).idx']);
    t_str_match(r.cpf.events(1).name, 'VLIM', [t 'events(1).name']);
    
    t = 'CPF to nose pt (pseudo arc length) w/violated V lims: ';
    mpopt_vlim = mpoption(mpopt, 'cpf.stop_at', 'NOSE', 'cpf.parameterization', 3,'cpf.enforce_v_lims',1);
    mpopt_vlim = mpoption(mpopt_vlim, 'cpf.adapt_step', 1);
    mpcb1 = mpcb;
    mpcb1.bus(6, VMIN) = 0.98;
    r = runcpf(mpcb1, mpct, mpopt_vlim);
    iterations = 1;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0, 6, [t 'max_lam']);
    t_is(r.gen(pvg, PG), Pg_lam(r.cpf.max_lam), 12, [t, 'Pg']);
    t_is(r.bus(:, PD), Pd_lam(r.cpf.max_lam), 12, [t, 'Pd']);
    t_is(size(r.cpf.V_hat), [10 iterations], 12, [t 'size(V_hat)']);
    t_is(size(r.cpf.V), [10 iterations], 12, [t 'size(V)']);
    t_is(size(r.cpf.lam_hat), [1 iterations], 12, [t 'size(lam_hat)']);
    t_is(size(r.cpf.lam), [1 iterations], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'voltage magnitude limit violated in base case'), [t 'done_msg']);
    t_is(length(r.cpf.events), 0, 12, [t 'length(events) == 0']);
    
    t = 'CPF to nose pt (pseudo arc length) w/V+flow lims: ';
    mpopt_vlim = mpoption(mpopt, 'cpf.stop_at', 'NOSE', 'cpf.parameterization', 3,'cpf.enforce_v_lims',1,'cpf.enforce_flow_lims',1);
    mpopt_vlim = mpoption(mpopt_vlim, 'cpf.adapt_step', 1);
    r = runcpf(mpcb, mpct, mpopt_vlim);
    iterations = 3;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0.110684, 6, [t 'max_lam']);
    t_is(r.gen(pvg, PG), Pg_lam(r.cpf.max_lam), 12, [t, 'Pg']);
    t_is(r.bus(:, PD), Pd_lam(r.cpf.max_lam), 12, [t, 'Pd']);
    t_is(size(r.cpf.V_hat), [10 iterations+1], 12, [t 'size(V_hat)']);
    t_is(size(r.cpf.V), [10 iterations+1], 12, [t 'size(V)']);
    t_is(size(r.cpf.lam_hat), [1 iterations+1], 12, [t 'size(lam_hat)']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'branch flow limit reached'), [t 'done_msg']);
    t_is(length(r.cpf.events), 1, 12, [t 'length(events) == 1']);
    t_is(r.cpf.events(1).k, iterations, 12, [t 'events(1).k']);
    t_is(r.cpf.events(1).idx, 5, 12, [t 'events(1).idx']);
    t_str_match(r.cpf.events(1).name, 'FLIM', [t 'events(1).name']);
    
    t = 'CPF to nose pt (pseudo arc length) w/PQ lims: ';
    mpopt_pqlim = mpoption(mpopt, 'cpf.stop_at', 'NOSE', 'cpf.parameterization', 3,'cpf.enforce_q_lims',1,'cpf.enforce_p_lims',1);
    mpopt_pqlim = mpoption(mpopt_pqlim, 'cpf.adapt_step', 1);
    r = runcpf(mpcb, mpct, mpopt_pqlim);
    iterations = 20;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0.833343, 3, [t 'max_lam']);
    ePg = Pg_lam(r.cpf.max_lam);
    t_is(r.gen(pvg(1), PG), ePg(1), 12, [t, 'Pg']);
    t_is(r.bus(:, PD), Pd_lam(r.cpf.max_lam), 12, [t, 'Pd']);
    t_is(size(r.cpf.V_hat), [10 iterations+1], 12, [t 'size(V_hat)']);
    t_is(size(r.cpf.V), [10 iterations+1], 12, [t 'size(V)']);
    t_is(size(r.cpf.lam_hat), [1 iterations+1], 12, [t 'size(lam_hat)']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'Reached steady state loading limit'), [t 'done_msg']);
    t_is(r.bus(:, BUS_TYPE), [1;3;2;4;1;1;1;1;1;1], 12, [t 'bus type']);
    ek = [7 14 iterations];
    eidx = [3 1 1];
    ename = {'PLIM', 'QLIM', 'NOSE'};
    ne = length(ek);
    t_is(length(r.cpf.events), ne, 12, sprintf('%ssize(events) == %d', t, ne));
    for j = 1:ne
        t_is(r.cpf.events(j).k, ek(j), 12, sprintf('%sevents(%d).k == %d', t, j, ek(j)));
        t_is(r.cpf.events(j).idx, eidx(j), 12, sprintf('%sevents(%d).idx == %d', t, j, eidx(j)));
        t_str_match(r.cpf.events(j).name, ename{j}, sprintf('%sevents(%d).name = ''%s''', t, j, ename{j}));
    end

    t = 'CPF (full trace) (arc length) : ';
    mpopt = mpoption(mpopt, 'cpf.stop_at', 'FULL', 'cpf.parameterization', 2);
    r = runcpf(mpcb, mpct, mpopt);
    iterations = 48;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0.99025, 3, [t 'max_lam']);
    t_is(r.gen(pvg, PG), Pg_lam(0), 12, [t, 'Pg']);
    t_is(r.bus(:, PD), Pd_lam(0), 12, [t, 'Pd']);
    t_is(size(r.cpf.V_hat), [10 iterations+1], 12, [t 'size(V_hat)']);
    t_is(size(r.cpf.V), [10 iterations+1], 12, [t 'size(V)']);
    t_is(size(r.cpf.lam_hat), [1 iterations+1], 12, [t 'size(lam_hat)']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'Traced full continuation curve in'), [t 'done_msg']);
    t_is(length(r.cpf.events), 1, 12, [t 'length(events) == 1']);
    t_is(r.cpf.events(1).k, iterations, 12, [t 'events(1).k']);
    t_is(r.cpf.events(1).idx, 1, 12, [t 'events(1).idx']);
    t_str_match(r.cpf.events(1).name, 'TARGET_LAM', [t 'events(1).name']);

    t = 'CPF (full trace) (pseudo arc length) : ';
    mpopt = mpoption(mpopt, 'cpf.stop_at', 'FULL', 'cpf.parameterization', 3);
    r = runcpf(mpcb, mpct, mpopt);
    iterations = 48;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0.99025, 3, [t 'max_lam']);
    t_is(r.gen(pvg, PG), Pg_lam(0), 12, [t, 'Pg']);
    t_is(r.bus(:, PD), Pd_lam(0), 12, [t, 'Pd']);
    t_is(size(r.cpf.V_hat), [10 iterations+1], 12, [t 'size(V_hat)']);
    t_is(size(r.cpf.V), [10 iterations+1], 12, [t 'size(V)']);
    t_is(size(r.cpf.lam_hat), [1 iterations+1], 12, [t 'size(lam_hat)']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'Traced full continuation curve in'), [t 'done_msg']);
    t_is(length(r.cpf.events), 1, 12, [t 'length(events) == 1']);
    t_is(r.cpf.events(1).k, iterations, 12, [t 'events(1).k']);
    t_is(r.cpf.events(1).idx, 1, 12, [t 'events(1).idx']);
    t_str_match(r.cpf.events(1).name, 'TARGET_LAM', [t 'events(1).name']);

    t = 'CPF (full trace) (pseudo arc length) w/Q lims: ';
    mpopt_qlim = mpoption(mpopt, 'cpf.stop_at', 'FULL', 'cpf.parameterization', 3,'cpf.enforce_q_lims',1);
    mpopt_qlim = mpoption(mpopt_qlim, 'cpf.adapt_step', 1);
    r = runcpf(mpcb, mpct, mpopt_qlim);
    iterations = 43;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0.795628, 6, [t 'max_lam']);
    ePg = Pg_lam(r.cpf.lam(end));
    t_is(r.gen(:, PG), [178.8688949; -64.7387022; 366.2496127], 6, [t, 'Pg']);
%r.gen(pvg, PG), Pg_lam(r.cpf.lam(end))
    t_is(r.bus(:, PD), Pd_lam(r.cpf.lam(end)), 12, [t, 'Pd']);
    t_is(size(r.cpf.V_hat), [10 iterations+1], 12, [t 'size(V_hat)']);
    t_is(size(r.cpf.V), [10 iterations+1], 12, [t 'size(V)']);
    t_is(size(r.cpf.lam_hat), [1 iterations+1], 12, [t 'size(lam_hat)']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_ok(regexp(r.cpf.done_msg, 'No REF or PV (buses|nodes) remaining.'), [t 'done_msg']);
    t_is(r.bus(:, BUS_TYPE), [1;1;3;4;1;1;1;1;1;1], 12, [t 'bus type']);
    ek = [12 22 iterations];
    eidx = [1 3 2];
    ename = {'QLIM', 'QLIM', 'QLIM'};
    ne = length(ek);
    t_is(length(r.cpf.events), ne, 12, sprintf('%ssize(events) == %d', t, ne));
    for j = 1:ne
        t_is(r.cpf.events(j).k, ek(j), 12, sprintf('%sevents(%d).k == %d', t, j, ek(j)));
        t_is(r.cpf.events(j).idx, eidx(j), 12, sprintf('%sevents(%d).idx == %d', t, j, eidx(j)));
        t_str_match(r.cpf.events(j).name, ename{j}, sprintf('%sevents(%d).name = ''%s''', t, j, ename{j}));
    end

    t = 'bug #12 : early termination : ';
    mpcbx = mpcb;
    mpctx = mpct;
    mpcbx.gen(1, QMAX) = 24.07;
    mpctx.gen(1, QMAX) = 24.07;
    r = runcpf(mpcbx, mpctx, mpopt_qlim);
    iterations = 40;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);

    t = 'all buses isolated : ';
    mpcbx.bus(:, BUS_TYPE) = NONE;
    try
        r = runcpf(mpcbx, mpctx, mpopt);
        t_is(r.success, 0, 12, [t 'success = 0']);
    catch
        t_ok(0, [t 'unexpected fatal error']);
    end

    t = 'issues #23,25 : direction switching : ';
    if have_feature('octave') && have_feature('octave', 'vnum') < 4.004
        %% eigs(-,-,'sr',-) throws error if all eigenvalues are real
        t_skip(10, 'eigs issue w/Octave < 4.4')
    else
        mpcbx = loadcase('case4gs');
        mpcbx.bus(4, [PD QD]) = [-1820 -1307.2];
        mpcbx.gen(1, [PG QMAX QMIN PMAX]) = [300 200 -1e9 250];
        mpcbx.gen(2, [QMAX PMAX]) = [1e9 500];
        mpctx = mpcbx;
        mpctx.bus(3, [PD QD]) = 2 * mpcbx.bus(3, [PD QD]);
        r = runcpf(mpcbx, mpctx, mpopt_qlim);
        iterations = 46;
        t_ok(r.success, [t 'success']);
        t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
        t_is(r.cpf.max_lam, 2.98036879, 8, [t 'max_lam']);
        t_is(length(r.cpf.events), 2, 12, [t 'length(events) == 2']);
        t_is(r.cpf.events(1).k, 25, 12, [t 'events(1).k']);
        t_is(r.cpf.events(1).idx, 1, 12, [t 'events(1).idx']);
        t_str_match(r.cpf.events(1).name, 'QLIM', [t 'events(1).name']);
        t_is(r.cpf.events(2).k, iterations, 12, [t 'events(2).k']);
        t_is(r.cpf.events(2).idx, 1, 12, [t 'events(1=2).idx']);
        t_str_match(r.cpf.events(2).name, 'TARGET_LAM', [t 'events(2).name']);
    end

    if have_feature('mp_core')
        t_skip(88, 'user callbacks')
    else
    t = '1 user callback : ';
    mpopt = mpoption(mpopt, 'cpf.stop_at', 0.7, 'cpf.parameterization', 3);
    mpopt = mpoption(mpopt, 'cpf.adapt_step', 1);
    mpopt = mpoption(mpopt, 'cpf.user_callback', 't_cpf_cb1');
    r = runcpf(mpcb, mpct, mpopt);
    iterations = 9;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0.7, 12, [t 'max_lam']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'Reached desired lambda 0.7'), [t 'done_msg']);
    t_is(length(r.cpf.events), 1, 12, [t 'length(events) == 1']);
    t_is(r.cpf.events(1).k, iterations, 12, [t 'events(1).k']);
    t_is(r.cpf.events(1).idx, 1, 12, [t 'events(1).idx']);
    t_str_match(r.cpf.events(1).name, 'TARGET_LAM', [t 'events(1).name']);
    t_ok(isfield(r.cpf, 'cb1'), [t 'isfield cpf.cb1']);
    t_ok(isstruct(r.cpf.cb1), [t 'isstruct cpf.cb1']);
    t_ok(isfield(r.cpf.cb1, 'initial'), [t 'isfield cpf.cb1.initial']);
    t_ok(isfield(r.cpf.cb1, 'iteration'), [t 'isfield cpf.cb1.iteration']);
    t_ok(isfield(r.cpf.cb1, 'final'), [t 'isfield cpf.cb1.final']);
    t_is(r.cpf.cb1.initial, 1, 12, [t 'r.cpf.cb1.initial']);
    t_is(r.cpf.cb1.iteration, iterations, 12, [t 'r.cpf.cb1.iterations']);
    t_is(r.cpf.cb1.final, 1, 12, [t 'r.cpf.cb1.final']);
    t_str_match(r.cpf.shared, '1111111111', [t 'r.cpf.shared']);

    t = '1 user callback : ';
    cb1 = struct('fcn', 't_cpf_cb1', 'priority', 10);
    mpopt = mpoption(mpopt, 'cpf.stop_at', 0.7, 'cpf.parameterization', 3);
    mpopt = mpoption(mpopt, 'cpf.adapt_step', 1);
    mpopt = mpoption(mpopt, 'cpf.user_callback', cb1);
    r = runcpf(mpcb, mpct, mpopt);
    iterations = 9;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0.7, 12, [t 'max_lam']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'Reached desired lambda 0.7'), [t 'done_msg']);
    t_is(length(r.cpf.events), 1, 12, [t 'length(events) == 1']);
    t_is(r.cpf.events(1).k, iterations, 12, [t 'events(1).k']);
    t_is(r.cpf.events(1).idx, 1, 12, [t 'events(1).idx']);
    t_str_match(r.cpf.events(1).name, 'TARGET_LAM', [t 'events(1).name']);
    t_ok(isfield(r.cpf, 'cb1'), [t 'isfield cpf.cb1']);
    t_ok(isstruct(r.cpf.cb1), [t 'isstruct cpf.cb1']);
    t_ok(isfield(r.cpf.cb1, 'initial'), [t 'isfield cpf.cb1.initial']);
    t_ok(isfield(r.cpf.cb1, 'iteration'), [t 'isfield cpf.cb1.iteration']);
    t_ok(isfield(r.cpf.cb1, 'final'), [t 'isfield cpf.cb1.final']);
    t_is(r.cpf.cb1.initial, 1, 12, [t 'r.cpf.cb1.initial']);
    t_is(r.cpf.cb1.iteration, iterations, 12, [t 'r.cpf.cb1.iterations']);
    t_is(r.cpf.cb1.final, 1, 12, [t 'r.cpf.cb1.final']);
    t_str_match(r.cpf.shared, '1111111111', [t 'r.cpf.shared']);

    t = '2 user callbacks (with args) : ';
    cb_args = struct('initial', 20, 'iteration', 2, 'final', 200);
    cb2 = struct('fcn', 't_cpf_cb2', 'args', cb_args);
    mpopt = mpoption(mpopt, 'cpf.user_callback', {'t_cpf_cb1', cb2});
    r = runcpf(mpcb, mpct, mpopt);
    iterations = 9;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0.7, 12, [t 'max_lam']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'Reached desired lambda 0.7'), [t 'done_msg']);
    t_is(length(r.cpf.events), 1, 12, [t 'length(events) == 1']);
    t_is(r.cpf.events(1).k, iterations, 12, [t 'events(1).k']);
    t_is(r.cpf.events(1).idx, 1, 12, [t 'events(1).idx']);
    t_str_match(r.cpf.events(1).name, 'TARGET_LAM', [t 'events(1).name']);
    t_ok(isfield(r.cpf, 'cb1'), [t 'isfield cpf.cb1']);
    t_ok(isstruct(r.cpf.cb1), [t 'isstruct cpf.cb1']);
    t_ok(isfield(r.cpf.cb1, 'initial'), [t 'isfield cpf.cb1.initial']);
    t_ok(isfield(r.cpf.cb1, 'iteration'), [t 'isfield cpf.cb1.iteration']);
    t_ok(isfield(r.cpf.cb1, 'final'), [t 'isfield cpf.cb1.final']);
    t_is(r.cpf.cb1.initial, 1, 12, [t 'r.cpf.cb1.initial']);
    t_is(r.cpf.cb1.iteration, iterations, 12, [t 'r.cpf.cb1.iterations']);
    t_is(r.cpf.cb1.final, 1, 12, [t 'r.cpf.cb1.final']);
    t_ok(isfield(r.cpf, 'cb2'), [t 'isfield cpf.cb2']);
    t_ok(isstruct(r.cpf.cb2), [t 'isstruct cpf.cb2']);
    t_ok(isfield(r.cpf.cb2, 'initial'), [t 'isfield cpf.cb2.initial']);
    t_ok(isfield(r.cpf.cb2, 'iteration'), [t 'isfield cpf.cb2.iteration']);
    t_ok(isfield(r.cpf.cb2, 'final'), [t 'isfield cpf.cb2.final']);
    t_is(r.cpf.cb2.initial, 20, 12, [t 'r.cpf.cb2.initial']);
    t_is(r.cpf.cb2.iteration, 2*iterations, 12, [t 'r.cpf.cb2.iterations']);
    t_is(r.cpf.cb2.final, 200, 12, [t 'r.cpf.cb2.final']);
    t_str_match(r.cpf.shared, '12121212121212121212', [t 'r.cpf.shared']);

    t = '2 user callbacks (with priority & args) : ';
    cb_args = struct('initial', 20, 'iteration', 2, 'final', 200);
    cb2 = struct('fcn', 't_cpf_cb2', 'priority', 21, 'args', cb_args);
    mpopt = mpoption(mpopt, 'cpf.user_callback', {'t_cpf_cb1', cb2});
    r = runcpf(mpcb, mpct, mpopt);
    iterations = 9;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0.7, 12, [t 'max_lam']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'Reached desired lambda 0.7'), [t 'done_msg']);
    t_is(length(r.cpf.events), 1, 12, [t 'length(events) == 1']);
    t_is(r.cpf.events(1).k, iterations, 12, [t 'events(1).k']);
    t_is(r.cpf.events(1).idx, 1, 12, [t 'events(1).idx']);
    t_str_match(r.cpf.events(1).name, 'TARGET_LAM', [t 'events(1).name']);
    t_ok(isfield(r.cpf, 'cb1'), [t 'isfield cpf.cb1']);
    t_ok(isstruct(r.cpf.cb1), [t 'isstruct cpf.cb1']);
    t_ok(isfield(r.cpf.cb1, 'initial'), [t 'isfield cpf.cb1.initial']);
    t_ok(isfield(r.cpf.cb1, 'iteration'), [t 'isfield cpf.cb1.iteration']);
    t_ok(isfield(r.cpf.cb1, 'final'), [t 'isfield cpf.cb1.final']);
    t_is(r.cpf.cb1.initial, 1, 12, [t 'r.cpf.cb1.initial']);
    t_is(r.cpf.cb1.iteration, iterations, 12, [t 'r.cpf.cb1.iterations']);
    t_is(r.cpf.cb1.final, 1, 12, [t 'r.cpf.cb1.final']);
    t_ok(isfield(r.cpf, 'cb2'), [t 'isfield cpf.cb2']);
    t_ok(isstruct(r.cpf.cb2), [t 'isstruct cpf.cb2']);
    t_ok(isfield(r.cpf.cb2, 'initial'), [t 'isfield cpf.cb2.initial']);
    t_ok(isfield(r.cpf.cb2, 'iteration'), [t 'isfield cpf.cb2.iteration']);
    t_ok(isfield(r.cpf.cb2, 'final'), [t 'isfield cpf.cb2.final']);
    t_is(r.cpf.cb2.initial, 20, 12, [t 'r.cpf.cb2.initial']);
    t_is(r.cpf.cb2.iteration, 2*iterations, 12, [t 'r.cpf.cb2.iterations']);
    t_is(r.cpf.cb2.final, 200, 12, [t 'r.cpf.cb2.final']);
    t_str_match(r.cpf.shared, '21212121212121212121', [t 'r.cpf.shared']);
    end  %% if mp_core

    t = 'case300 w/Q lims : ';
    mpopt = mpoption('out.all', 0, 'verbose', verbose);
    mpopt = mpoption(mpopt, 'cpf.stop_at', 'FULL', 'cpf.step', 0.1);
    mpopt = mpoption(mpopt, 'cpf.step_max', 0.2);
    mpopt = mpoption(mpopt, 'cpf.adapt_step', 0);
    mpopt = mpoption(mpopt, 'cpf.enforce_q_lims', 1);
    mpcb = loadcase('case300');                         % load base case
    mpct = mpcb;                                        % set up target case with
    mpct.gen(:, [PG QG]) = mpcb.gen(:, [PG QG]) * 2.5;  % increased generation
    mpct.bus(:, [PD QD]) = mpcb.bus(:, [PD QD]) * 2.5;  % and increased load

    r = runcpf(mpcb, mpct, mpopt);
    iterations = 43;
    t_ok(r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 0.038307, 5, [t 'max_lam']);
    t_is(size(r.cpf.V_hat), [300 iterations+1], 12, [t 'size(V_hat)']);
    t_is(size(r.cpf.V), [300 iterations+1], 12, [t 'size(V)']);
    t_is(size(r.cpf.lam_hat), [1 iterations+1], 12, [t 'size(lam_hat)']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'Traced full continuation curve in'), [t 'done_msg']);
    ebt = ones(size(r.bus(:, BUS_TYPE)));
    ebt([69 76 77 80 88 98 117 122 132 155 164:166 169 170 177 192 199:201 206 209 212 217 218 220:222 247 248 250 253 254 263:265]) = 2;
    ebt(69) = 3;
    t_is(r.bus(:, BUS_TYPE), ebt, 12, [t 'bus type']);
    qmin_violation = max(0, r.gen(:, QMIN) - r.gen(:, QG));
    qmax_violation = max(0, r.gen(:, QG) - r.gen(:, QMAX));
    t_is(qmin_violation, 0, 12, [t 'Qmin violation']);
    t_is(qmax_violation, 0, 12, [t 'Qmax violation']);
    ek = [1 3 5 12 14 17 19 22 24 28 38 41 43];
    eidx = [26 19 67 15 59 69 50 66 54 68 87 89 1];
    ename = {'QLIM', 'QLIM', 'QLIM', 'QLIM', 'QLIM', 'QLIM', 'QLIM', 'QLIM', 'QLIM', 'QLIM', 'QLIM', 'QLIM', 'TARGET_LAM'};
    ne = length(ek);
    t_is(length(r.cpf.events), ne, 12, sprintf('%ssize(events) == %d', t, ne));
    for j = 1:ne
        t_is(r.cpf.events(j).k, ek(j), 12, sprintf('%sevents(%d).k == %d', t, j, ek(j)));
        t_is(r.cpf.events(j).idx, eidx(j), 12, sprintf('%sevents(%d).idx == %d', t, j, eidx(j)));
        t_str_match(r.cpf.events(j).name, ename{j}, sprintf('%sevents(%d).name = ''%s''', t, j, ename{j}));
    end

    t = 'case14 (unsuccessful) : ';
    if have_feature('mp_core')
        t_skip(9, [t 'MP-Core not affected by angle wrap-around failure']);
    else
    mpopt = mpoption(mpopt, 'cpf.adapt_step', 1);
    mpopt = mpoption(mpopt, 'cpf.enforce_q_lims', 0);
    mpcb = loadcase('case14');                          % load base case
    mpct = mpcb;                                        % set up target case with
    mpct.gen(:, [PG QG]) = mpcb.gen(:, [PG QG]) * 2.5;  % increased generation
    mpct.bus(:, [PD QD]) = mpcb.bus(:, [PD QD]) * 2.5;  % and increased load

    r = runcpf(mpcb, mpct, mpopt);
    iterations = 79;
    t_ok(~r.success, [t 'success']);
    t_is(r.cpf.iterations, iterations, 12, [t 'iterations']);
    t_is(r.cpf.max_lam, 2.0401678, 6, [t 'max_lam']);
    t_is(size(r.cpf.V_hat), [14 iterations+1], 12, [t 'size(V_hat)']);
    t_is(size(r.cpf.V), [14 iterations+1], 12, [t 'size(V)']);
    t_is(size(r.cpf.lam_hat), [1 iterations+1], 12, [t 'size(lam_hat)']);
    t_is(size(r.cpf.lam), [1 iterations+1], 12, [t 'size(lam)']);
    t_ok(strfind(r.cpf.done_msg, 'Corrector did not converge in'), [t 'done_msg']);
    t_is(length(r.cpf.events), 0, 12, [t 'length(events) == 0']);
    end

    if have_feature('octave')
        warning(s1.state, file_in_path_warn_id);
    end
end

t_end;

if nargout
    res = r;
end
