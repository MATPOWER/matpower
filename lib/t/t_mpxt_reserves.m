function t_mpxt_reserves(quiet)
% t_mpxt_reserves - Tests mp.xt_reserves extension.

%   MATPOWER
%   Copyright (c) 2009-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

num_tests = 66;
diff_tool = 'bbdiff';
show_diff_on_fail = false;

t_begin(num_tests, quiet);

if quiet
    verbose = 0;
else
    verbose = 0;
end
if have_feature('octave')
    file_in_path_warn_id = 'Octave:data-file-in-path';
    s1 = warning('query', file_in_path_warn_id);
    warning('off', file_in_path_warn_id);
    near_sing_matrix_warn_id = 'Octave:nearly-singular-matrix';
else
    near_sing_matrix_warn_id = 'MATLAB:nearlySingularMatrix';
end
s = warning('query', near_sing_matrix_warn_id);
warning('off', near_sing_matrix_warn_id);

[pathstr, name, ext] = fileparts(which('t_pretty_print'));

casefile = 't_case30_userfcns';
mpopt = mpoption('opf.violation', 1e-6, 'mips.gradtol', 1e-8, ...
        'mips.comptol', 1e-8, 'mips.costtol', 1e-9);
mpopt = mpoption(mpopt, 'out.all', 0, 'verbose', verbose, 'opf.ac.solver', 'MIPS');
mpopt = mpoption(mpopt, 'opf.dc.solver', 'MIPS');
mpopt = mpoption(mpopt, 'out.gen', 1);
mpopt = mpoption(mpopt, 'out.lim.all', 2);
% mpopt = mpoption(mpopt, 'out.lim.all', -1, 'out.sys_sum', 0, 'out.bus', 0, 'out.branch', 0, 'out.gen', 1, 'out.lim.v', 0, 'out.lim.pg', 0, 'out.lim.qg', 0, 'out.lim.line', 0); mpopt.out.load = 0; mpopt.out.shunt = 0; mpopt.out.lim.elm.reserve_gen = 2;

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

t = 'run_opf(''t_case30_userfcns'', ...) : ';
rn = fix(1e9*rand);
fname = sprintf('pp_mpxt_reserve_%d', 1);
fname_e = fullfile(pathstr, 'pretty-printing', sprintf('%s.txt', fname));
fname_g = sprintf('%s_%d.txt', fname, rn);
task = run_opf(casefile, mpopt, 'mpx', mp.xt_reserves, 'print_fname', fname_g);
rg = task.dm.elements.reserve_gen.tab;
rz = task.dm.elements.reserve_zone.tab;
t_ok(task.success, [t 'success']);
t_is(rg.r, [25; 15; 0; 0; 19.3906; 0.6094], 4, [t 'R']);
t_is(rg.prc, [2; 2; 2; 2; 5.5; 5.5], 6, [t 'prc']);
t_is(rg.mu_lb, [0; 0; 1; 2; 0; 0], 7, [t 'mu.l']);
t_is(rg.mu_ub, [0.1; 0; 0; 0; 0; 0], 7, [t 'mu.u']);
t_is(rg.mu_pg_ub, [0; 0; 0; 0; 0.5; 0], 7, [t 'mu.Pmax']);
mpc = loadcase(casefile);
t_is(rg.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(rg.qty, mpc.reserves.qty, 12, [t 'qty']);
t_is(sum(rg.total_cost), 177.8047, 4, [t 'totalcost']);
reps = { {' -(0.0+) ', '  $1 ', 1, 1}, ...
         {'in (.*) seconds \((.*) setup \+ (.*) solve\)', ...
            'in 0.00 seconds (0.00 setup + 0.00 solve)', 1, 1} };
if ~t_file_match(fname_g, fname_e, [t 'pretty printing'], reps, 1);
    fprintf('  compare these 2 files:\n    %s\n    %s\n', fname_g, fname_e);
    if show_diff_on_fail
        cmd = sprintf('%s %s %s', diff_tool, fname_g, fname_e);
        [status, result] = system(cmd);
        keyboard
    end
end

t = 'gen 5 no reserves : ';
mpc = loadcase(casefile);
mpc.reserves.zones(:, 5) = 0;
mpc.reserves.cost(5) = [];
mpc.reserves.qty(5) = [];
rn = fix(1e9*rand);
fname = sprintf('pp_mpxt_reserve_%d', 2);
fname_e = fullfile(pathstr, 'pretty-printing', sprintf('%s.txt', fname));
fname_g = sprintf('%s_%d.txt', fname, rn);
task = run_opf(mpc, mpopt, 'mpx', mp.xt_reserves, 'print_fname', fname_g);
rg = task.dm.elements.reserve_gen.tab;
rz = task.dm.elements.reserve_zone.tab;
t_ok(task.success, [t 'success']);
t_is(rg.r, [25; 15; 0; 0; 20], 4, [t 'R']);
t_is(rg.prc, [2; 2; 2; 2; 5.5], 6, [t 'prc']);
t_is(rg.mu_lb, [0; 0; 1; 2; 0], 7, [t 'mu.l']);
t_is(rg.mu_ub, [0.1; 0; 0; 0; 0], 6, [t 'mu.u']);
t_is(rg.mu_pg_ub, [0; 0; 0; 0; 0], 7, [t 'mu.Pmax']);
t_is(rg.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(rg.qty, mpc.reserves.qty, 12, [t 'qty']);
t_is(sum(rg.total_cost), 187.5, 4, [t 'totalcost']);
reps = { {' -(0.0+) ', '  $1 ', 1, 1}, ...
         {'in (.*) seconds \((.*) setup \+ (.*) solve\)', ...
            'in 0.00 seconds (0.00 setup + 0.00 solve)', 1, 1} };
if ~t_file_match(fname_g, fname_e, [t 'pretty printing'], reps, 1);
    fprintf('  compare these 2 files:\n    %s\n    %s\n', fname_g, fname_e);
    if show_diff_on_fail
        cmd = sprintf('%s %s %s', diff_tool, fname_g, fname_e);
        [status, result] = system(cmd);
        keyboard
    end
end

t = 'extra offline gen : ';
mpc = loadcase(casefile);
idx = [1:3 5 4:6];
mpc.gen = mpc.gen(idx, :);
mpc.gencost = mpc.gencost(idx, :);
mpc.reserves.zones = mpc.reserves.zones(:, idx);
mpc.reserves.cost = mpc.reserves.cost(idx);
mpc.reserves.qty = mpc.reserves.qty(idx);
mpc.gen(4, GEN_STATUS) = 0;
task = run_opf(mpc, mpopt, 'mpx', mp.xt_reserves);
rg = task.dm.elements.reserve_gen.tab;
rz = task.dm.elements.reserve_zone.tab;
t_ok(task.success, [t 'success']);
t_is(rg.r, [25; 15; 0; 0; 0; 19.3906; 0.6094], 4, [t 'R']);
t_is(rg.prc, [2; 2; 2; 5.5; 2; 5.5; 5.5], 6, [t 'prc']);
t_is(rg.mu_lb, [0; 0; 1; 0; 2; 0; 0], 7, [t 'mu.l']);
t_is(rg.mu_ub, [0.1; 0; 0; 0; 0; 0; 0], 7, [t 'mu.u']);
t_is(rg.mu_pg_ub, [0; 0; 0; 0; 0; 0.5; 0], 7, [t 'mu.Pmax']);
t_is(rg.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(rg.qty, mpc.reserves.qty, 12, [t 'qty']);
t_is(sum(rg.total_cost), 177.8047, 4, [t 'totalcost']);

t = 'both extra & gen 6 no res : ';
mpc = loadcase(casefile);
idx = [1:3 5 4:6];
mpc.gen = mpc.gen(idx, :);
mpc.gencost = mpc.gencost(idx, :);
mpc.reserves.zones = mpc.reserves.zones(:, idx);
mpc.reserves.cost = mpc.reserves.cost(idx);
mpc.reserves.qty = mpc.reserves.qty(idx);
mpc.gen(4, GEN_STATUS) = 0;
mpc.reserves.zones(:, 6) = 0;
mpc.reserves.cost(6) = [];
mpc.reserves.qty(6) = [];
rn = fix(1e9*rand);
fname = sprintf('pp_mpxt_reserve_%d', 3);
fname_e = fullfile(pathstr, 'pretty-printing', sprintf('%s.txt', fname));
fname_g = sprintf('%s_%d.txt', fname, rn);
task = run_opf(mpc, mpopt, 'mpx', mp.xt_reserves, 'print_fname', fname_g);
rg = task.dm.elements.reserve_gen.tab;
rz = task.dm.elements.reserve_zone.tab;
t_ok(task.success, [t 'success']);
t_is(rg.r, [25; 15; 0; 0; 0; 20], 4, [t 'R']);
t_is(rg.prc, [2; 2; 2; 5.5; 2; 5.5], 6, [t 'prc']);
t_is(rg.mu_lb, [0; 0; 1; 0; 2; 0], 7, [t 'mu.l']);
t_is(rg.mu_ub, [0.1; 0; 0; 0; 0; 0], 6, [t 'mu.u']);
t_is(rg.mu_pg_ub, [0; 0; 0; 0; 0; 0], 7, [t 'mu.Pmax']);
t_is(rg.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(rg.qty, mpc.reserves.qty, 12, [t 'qty']);
t_is(sum(rg.total_cost), 187.5, 4, [t 'totalcost']);
reps = { {' -(0.0+) ', '  $1 ', 1, 1}, ...
         {'in (.*) seconds \((.*) setup \+ (.*) solve\)', ...
            'in 0.00 seconds (0.00 setup + 0.00 solve)', 1, 1} };
if ~t_file_match(fname_g, fname_e, [t 'pretty printing'], reps, 1);
    fprintf('  compare these 2 files:\n    %s\n    %s\n', fname_g, fname_e);
    if show_diff_on_fail
        cmd = sprintf('%s %s %s', diff_tool, fname_g, fname_e);
        [status, result] = system(cmd);
        keyboard
    end
end

t = 'no qty (Rmax) : ';
mpc = loadcase(casefile);
mpc.reserves = rmfield(mpc.reserves, 'qty');
task = run_opf(mpc, mpopt, 'mpx', mp.xt_reserves);
rg = task.dm.elements.reserve_gen.tab;
rz = task.dm.elements.reserve_zone.tab;
t_ok(task.success, [t 'success']);
t_is(rg.r, [39.3826; 0.6174; 0; 0; 19.3818; 0.6182], 4, [t 'R']);
t_is(rg.prc, [2; 2; 2; 2; 5.5; 5.5], 5, [t 'prc']);
t_is(rg.mu_lb, [0; 0; 1; 2; 0; 0], 5, [t 'mu.l']);
t_is(rg.mu_ub, [0; 0; 0; 0; 0; 0], 7, [t 'mu.u']);
t_is(rg.mu_pg_ub, [0.1; 0; 0; 0; 0.5; 0], 5, [t 'mu.Pmax']);
t_is(rg.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(sum(rg.total_cost), 176.3708, 4, [t 'totalcost']);

t = 'RAMP_10, no qty (Rmax) : ';
mpc = loadcase(casefile);
mpc.reserves = rmfield(mpc.reserves, 'qty');
mpc.gen(1, RAMP_10) = 25;
rn = fix(1e9*rand);
fname = sprintf('pp_mpxt_reserve_%d', 4);
fname_e = fullfile(pathstr, 'pretty-printing', sprintf('%s.txt', fname));
fname_g = sprintf('%s_%d.txt', fname, rn);
task = run_opf(mpc, mpopt, 'mpx', mp.xt_reserves, 'print_fname', fname_g);
rg = task.dm.elements.reserve_gen.tab;
rz = task.dm.elements.reserve_zone.tab;
t_ok(task.success, [t 'success']);
t_is(rg.r, [25; 15; 0; 0; 19.3906; 0.6094], 4, [t 'R']);
t_is(rg.prc, [2; 2; 2; 2; 5.5; 5.5], 6, [t 'prc']);
t_is(rg.mu_lb, [0; 0; 1; 2; 0; 0], 7, [t 'mu.l']);
t_is(rg.mu_ub, [0.1; 0; 0; 0; 0; 0], 7, [t 'mu.u']);
t_is(rg.mu_pg_ub, [0; 0; 0; 0; 0.5; 0], 7, [t 'mu.Pmax']);
t_is(rg.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(sum(rg.total_cost), 177.8047, 4, [t 'totalcost']);
reps = { {' -(0.0+) ', '  $1 ', 1, 1}, ...
         {'in (.*) seconds \((.*) setup \+ (.*) solve\)', ...
            'in 0.00 seconds (0.00 setup + 0.00 solve)', 1, 1} };
if ~t_file_match(fname_g, fname_e, [t 'pretty printing'], reps, 1);
    fprintf('  compare these 2 files:\n    %s\n    %s\n', fname_g, fname_e);
    if show_diff_on_fail
        cmd = sprintf('%s %s %s', diff_tool, fname_g, fname_e);
        [status, result] = system(cmd);
        keyboard
    end
end

t = 'DC OPF : ';
mpc = loadcase(casefile);
mpopt = mpoption(mpopt, 'model', 'DC', 'out.sys_sum', 0);
rn = fix(1e9*rand);
fname = sprintf('pp_mpxt_reserve_%d', 5);
fname_e = fullfile(pathstr, 'pretty-printing', sprintf('%s.txt', fname));
fname_g = sprintf('%s_%d.txt', fname, rn);
task = run_opf(mpc, mpopt, 'mpx', mp.xt_reserves, 'print_fname', fname_g);
rg = task.dm.elements.reserve_gen.tab;
rz = task.dm.elements.reserve_zone.tab;
t_ok(task.success, [t 'success']);
t_is(rg.r, [25; 15; 0; 0; 20; 0], 4, [t 'R']);
t_is(rg.prc, [2; 2; 2; 2; 5.33006533; 5.33006533], 6, [t 'prc']);
t_is(rg.mu_lb, [0; 0; 1; 2; 0; 0.169934666], 7, [t 'mu.l']);
t_is(rg.mu_ub, [0.1; 0; 0; 0; 0; 0], 7, [t 'mu.u']);
t_is(rg.mu_pg_ub, [0; 0; 0; 0; 0.33006533; 0], 7, [t 'mu.Pmax']);
t_is(rg.cost, mpc.reserves.cost, 12, [t 'cost']);
t_is(rg.qty, mpc.reserves.qty, 12, [t 'qty']);
t_is(sum(rg.total_cost), 177.5, 4, [t 'totalcost']);
reps = { {' -(0.0+) ', '  $1 ', 1, 1}, ...
         {'in (.*) seconds \((.*) setup \+ (.*) solve\)', ...
            'in 0.00 seconds (0.00 setup + 0.00 solve)', 1, 1} };
if ~t_file_match(fname_g, fname_e, [t 'pretty printing'], reps, 1);
    fprintf('  compare these 2 files:\n    %s\n    %s\n', fname_g, fname_e);
    if show_diff_on_fail
        cmd = sprintf('%s %s %s', diff_tool, fname_g, fname_e);
        [status, result] = system(cmd);
        keyboard
    end
end

if have_feature('octave')
    warning(s1.state, file_in_path_warn_id);
end
warning(s.state, near_sing_matrix_warn_id);

t_end;
