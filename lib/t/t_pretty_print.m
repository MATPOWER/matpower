function t_pretty_print(quiet)
% t_pretty_print - Tests for pretty printed output.

%   MATPOWER
%   Copyright (c) 2014-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

num_tests = 7;
diff_tool = 'bbdiff';
show_diff_on_fail = false;

t_begin(num_tests, quiet);

[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

if quiet
    verbose = 0;
else
    verbose = 0;
end
if have_feature('octave')
    file_in_path_warn_id = 'Octave:data-file-in-path';
    s1 = warning('query', file_in_path_warn_id);
    warning('off', file_in_path_warn_id);
end

[pathstr, name, ext] = fileparts(which('t_pretty_print'));

mpopt = mpoption('opf.violation', 1e-6, 'mips.gradtol', 1e-8, ...
        'mips.comptol', 1e-8, 'mips.costtol', 1e-9);
mpopt = mpoption(mpopt, 'out.all', 0, 'verbose', verbose);
mpopt = mpoption(mpopt, 'opf.ac.solver', 'MIPS');
mpopt = mpoption(mpopt, 'opf.dc.solver', 'MIPS');
% mpopt = mpoption(mpopt, 'out.sys_sum', 1);
% mpopt = mpoption(mpopt, 'out.area_sum', 1);
% mpopt = mpoption(mpopt, 'out.bus', 1);
% mpopt = mpoption(mpopt, 'out.branch', 1);
mpopt = mpoption(mpopt, 'out.gen', 1);
mpopt = mpoption(mpopt, 'out.lim.all', 2);

casefile = 'case9';

t = sprintf('run_pf(''%s'')', casefile);
rn = fix(1e9*rand);
fname = sprintf('pp_pf_%s', casefile);
fname_e = fullfile(pathstr, 'pretty-printing', sprintf('%s.txt', fname));
fname_g = sprintf('%s_%d.txt', fname, rn);
task = run_pf(casefile, mpopt, 'print_fname', fname_g);
reps = { {'in (.*) seconds \((.*) setup \+ (.*) solve\)', ...
            'in 0.00 seconds (0.00 setup + 0.00 solve)', 1, 1} };
if ~t_file_match(fname_g, fname_e, t, reps, 1);
    fprintf('  compare these 2 files:\n    %s\n    %s\n', fname_g, fname_e);
    if show_diff_on_fail
        cmd = sprintf('%s %s %s', diff_tool, fname_g, fname_e);
        [status, result] = system(cmd);
        keyboard
    end
end

t = sprintf('run_cpf({''%s'', ''%starget''})', casefile, casefile);
rn = fix(1e9*rand);
fname = sprintf('pp_cpf_%s', casefile);
fname_e = fullfile(pathstr, 'pretty-printing', sprintf('%s.txt', fname));
fname_g = sprintf('%s_%d.txt', fname, rn);
task = run_cpf({casefile, [casefile 'target']}, mpopt, 'print_fname', fname_g);
reps = { {'in (.*) seconds \((.*) setup \+ (.*) solve\)', ...
            'in 0.00 seconds (0.00 setup + 0.00 solve)', 1, 1} };
if ~t_file_match(fname_g, fname_e, t, reps, 1);
    fprintf('  compare these 2 files:\n    %s\n    %s\n', fname_g, fname_e);
    if show_diff_on_fail
        cmd = sprintf('%s %s %s', diff_tool, fname_g, fname_e);
        [status, result] = system(cmd);
        keyboard
    end
end

t = sprintf('run_opf(''%s'')', casefile);
rn = fix(1e9*rand);
fname = sprintf('pp_opf_%s', casefile);
fname_e = fullfile(pathstr, 'pretty-printing', sprintf('%s.txt', fname));
fname_g = sprintf('%s_%d.txt', fname, rn);
task = run_opf(casefile, mpopt, 'print_fname', fname_g);
reps = { {'in (.*) seconds \((.*) setup \+ (.*) solve\)', ...
            'in 0.00 seconds (0.00 setup + 0.00 solve)', 1, 1} };
if ~t_file_match(fname_g, fname_e, t, reps, 1);
    fprintf('  compare these 2 files:\n    %s\n    %s\n', fname_g, fname_e);
    if show_diff_on_fail
        cmd = sprintf('%s %s %s', diff_tool, fname_g, fname_e);
        [status, result] = system(cmd);
        keyboard
    end
end

%% load case
casefile = 't_auction_case';
mpc = loadcase(casefile);
mpc.gencost(1:5, 1:7) = ones(5, 1) * [2 0 0 3 0.1 45 0];
mpc.branch(18, RATE_A) = 25;

t = sprintf('run_pf(''%s'')', casefile);
rn = fix(1e9*rand);
fname = sprintf('pp_pf_%s', casefile);
fname_e = fullfile(pathstr, 'pretty-printing', sprintf('%s.txt', fname));
fname_g = sprintf('%s_%d.txt', fname, rn);
task = run_pf(mpc, mpopt, 'print_fname', fname_g);
reps = { {'in (.*) seconds \((.*) setup \+ (.*) solve\)', ...
            'in 0.00 seconds (0.00 setup + 0.00 solve)', 1, 1} };
if ~t_file_match(fname_g, fname_e, t, reps, 1);
    fprintf('  compare these 2 files:\n    %s\n    %s\n', fname_g, fname_e);
    if show_diff_on_fail
        cmd = sprintf('%s %s %s', diff_tool, fname_g, fname_e);
        [status, result] = system(cmd);
        keyboard
    end
end

t = sprintf('run_opf(''%s'')', casefile);
rn = fix(1e9*rand);
fname = sprintf('pp_opf_%s', casefile);
fname_e = fullfile(pathstr, 'pretty-printing', sprintf('%s.txt', fname));
fname_g = sprintf('%s_%d.txt', fname, rn);
task = run_opf(mpc, mpoption(mpopt, 'out.lim.all', 1), 'print_fname', fname_g);
reps = { {' -(0.0+) ', '  $1 ', 1, 1}, ...
         {'in (.*) seconds \((.*) setup \+ (.*) solve\)', ...
            'in 0.00 seconds (0.00 setup + 0.00 solve)', 1, 1} };
if ~t_file_match(fname_g, fname_e, t, reps, 1);
    fprintf('  compare these 2 files:\n    %s\n    %s\n', fname_g, fname_e);
    if show_diff_on_fail
        cmd = sprintf('%s %s %s', diff_tool, fname_g, fname_e);
        [status, result] = system(cmd);
        keyboard
    end
end

casefile = 't_case3p_a';

t = sprintf('run_pf(''%s'')', casefile);
rn = fix(1e9*rand);
fname = sprintf('pp_pf_%s', casefile);
fname_e = fullfile(pathstr, 'pretty-printing', sprintf('%s.txt', fname));
fname_g = sprintf('%s_%d.txt', fname, rn);
task = run_pf(casefile, mpopt, 'print_fname', fname_g, 'mpx', mp.xt_3p());
reps = { {'in (.*) seconds \((.*) setup \+ (.*) solve\)', ...
            'in 0.00 seconds (0.00 setup + 0.00 solve)', 1, 1} };
if ~t_file_match(fname_g, fname_e, t, reps, 1);
    fprintf('  compare these 2 files:\n    %s\n    %s\n', fname_g, fname_e);
    if show_diff_on_fail
        cmd = sprintf('%s %s %s', diff_tool, fname_g, fname_e);
        [status, result] = system(cmd);
        keyboard
    end
end

casefile = 't_case3p_g';

t = sprintf('run_opf(''%s'')', casefile);
rn = fix(1e9*rand);
fname = sprintf('pp_opf_%s', casefile);
fname_e = fullfile(pathstr, 'pretty-printing', sprintf('%s.txt', fname));
fname_g = sprintf('%s_%d.txt', fname, rn);
mpopt.exp.mpx = mp.xt_3p();
task = run_opf(casefile, mpopt, 'print_fname', fname_g);
reps = { {'in (.*) seconds \((.*) setup \+ (.*) solve\)', ...
            'in 0.00 seconds (0.00 setup + 0.00 solve)', 1, 1} };
if ~t_file_match(fname_g, fname_e, t, reps, 1);
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

t_end;
