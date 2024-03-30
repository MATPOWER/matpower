function t_printpf(quiet)
% t_printpf - Tests for printpf.

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

num_tests = 2;

t_begin(num_tests, quiet);

[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

casefile = 't_auction_case';
if quiet
    verbose = 0;
else
    verbose = 0;
end
if have_feature('octave')
    if have_feature('octave', 'vnum') >= 4
        file_in_path_warn_id = 'Octave:data-file-in-path';
    else
        file_in_path_warn_id = 'Octave:fopen-file-in-path';
    end
    s1 = warning('query', file_in_path_warn_id);
    warning('off', file_in_path_warn_id);
else
    sing_matrix_warn_id = 'MATLAB:nearlySingularMatrix';
    s2 = warning('query', sing_matrix_warn_id);
    warning('off', sing_matrix_warn_id);
end

t = 'printpf : ';
mpopt = mpoption('opf.violation', 1e-6, 'mips.gradtol', 1e-8, ...
        'mips.comptol', 1e-8, 'mips.costtol', 1e-9);
mpopt = mpoption(mpopt, 'out.all', 0, 'verbose', verbose);
mpopt = mpoption(mpopt, 'out.sys_sum', 1);
mpopt = mpoption(mpopt, 'out.area_sum', 1);
mpopt = mpoption(mpopt, 'out.bus', 1);
mpopt = mpoption(mpopt, 'out.branch', 1);
mpopt = mpoption(mpopt, 'out.gen', 1);
mpopt = mpoption(mpopt, 'out.lim.all', 2);
mpopt = mpoption(mpopt, 'opf.ac.solver', 'MIPS');
mpopt = mpoption(mpopt, 'opf.dc.solver', 'MIPS');

%% load case
mpc = loadcase(casefile);
mpc.gencost(1:5, 1:7) = ones(5, 1) * [2 0 0 3 0.1 45 0];
mpc.branch(18, RATE_A) = 25;

%% set up file names
fnameac = 'pretty_print_acopf.txt';
fnamedc = 'pretty_print_dcopf.txt';
fname = 'pretty_print';
rn = fix(1e9*rand);
tmpfnameac = sprintf('%s_acopf_%d.txt', fname, rn);
tmpfnamedc = sprintf('%s_dcopf_%d.txt', fname, rn);

r = runopf(mpc, mpopt, tmpfnameac);
reps = {{'Converged in (.*) seconds', 'Converged in 0.00 seconds'}, ...
        {' -0.0 ', '  0.0 ', 0}, ...
        {sprintf(' -0.0\n'), sprintf('  0.0\n'), 0}, ...
        {' -0.00 ', '  0.00 ', 0}, ...
        {sprintf(' -0.00\n'), sprintf('  0.00\n'), 0}, ...
        {' -0.000 ', '  0.000 ', 0}, ...
        {' -0.000*', '  0.000*', 0}};
t_file_match(tmpfnameac, fnameac, [t 'standard AC OPF'], reps, 1);

r = rundcopf(mpc, mpopt, tmpfnamedc);
reps = {{'Converged in (.*) seconds', 'Converged in 0.00 seconds'}, ...
        {' -0.0 ', '  0.0 ', 0}, ...
        {sprintf(' -0.0\n'), sprintf('  0.0\n'), 0}, ...
        {' -0.00 ', '  0.00 ', 0}, ...
        {sprintf(' -0.00\n'), sprintf('  0.00\n'), 0}, ...
        {' -0.000 ', '  0.000 ', 0}, ...
        {' -0.000*', '  0.000*', 0}, ...
        {'51.66 $/MWh @ bus 12', '51.66 $/MWh @ bus 13', 0}, ...
        {'51.66 $/MWh @ bus 16', '51.66 $/MWh @ bus 13', 0}, ...
        {'51.66 $/MWh @ bus 17', '51.66 $/MWh @ bus 13', 0}, ...
        {'53.05 $/MWh @ bus 18', '53.05 $/MWh @ bus 15', 0}, ...
        {'53.05 $/MWh @ bus 19', '53.05 $/MWh @ bus 15', 0}, ...
        {'53.05 $/MWh @ bus 20', '53.05 $/MWh @ bus 15', 0}};
t_file_match(tmpfnamedc, fnamedc, [t 'standard DC OPF'], reps, 1);

if have_feature('octave')
    warning(s1.state, file_in_path_warn_id);
else
    warning(s2.state, sing_matrix_warn_id);
end

t_end;
