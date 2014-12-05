function t_printpf(quiet)
%T_PRINTPF  Tests for PRINTPF.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2014 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

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
if have_fcn('octave')
    s1 = warning('query', 'Octave:load-file-in-path');
    warning('off', 'Octave:load-file-in-path');
end
s2 = warning('query', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:singularMatrix');

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
[pathstr, name, ext] = fileparts(which('t_printpf'));

fnameac = 'pretty_print_acopf.txt';
fnamedc = 'pretty_print_dcopf.txt';
fname = 'pretty_print';
rn = fix(1e9*rand);
tmpfnameac = sprintf('%s_acopf_%d.txt', fname, rn);
tmpfnamedc = sprintf('%s_dcopf_%d.txt', fname, rn);

r = runopf(mpc, mpopt, tmpfnameac);
got = fileread(tmpfnameac);
if size(got, 1) ~= 1    %% transpose if needed for Octave 3.4
    got = got';
end
got = regexprep(got, 'Converged in (.*) seconds', 'Converged in 0.00 seconds');
got = strrep(got, ' -0.00 ', '  0.00 ');
got = strrep(got, ' -0.000 ', '  0.000 ');
expected = fileread(fnameac);
if size(expected, 1) ~= 1   %% transpose if needed for Octave 3.4
    expected = expected';
end
t_ok(strcmp(got, expected), [t 'standard AC OPF']);
delete(tmpfnameac);

r = rundcopf(mpc, mpopt, tmpfnamedc);
got = fileread(tmpfnamedc);
if size(got, 1) ~= 1    %% transpose if needed for Octave 3.4
    got = got';
end
got = regexprep(got, 'Converged in (.*) seconds', 'Converged in 0.00 seconds');
got = strrep(got, ' -0.00 ', '  0.00 ');
got = strrep(got, ' -0.000 ', '  0.000 ');
got = strrep(got, '51.66 $/MWh @ bus 12', '51.66 $/MWh @ bus 13');
got = strrep(got, '51.66 $/MWh @ bus 16', '51.66 $/MWh @ bus 13');
got = strrep(got, '51.66 $/MWh @ bus 17', '51.66 $/MWh @ bus 13');
got = strrep(got, '53.05 $/MWh @ bus 18', '53.05 $/MWh @ bus 15');
got = strrep(got, '53.05 $/MWh @ bus 19', '53.05 $/MWh @ bus 15');
got = strrep(got, '53.05 $/MWh @ bus 20', '53.05 $/MWh @ bus 15');
expected = fileread(fnamedc);
if size(expected, 1) ~= 1   %% transpose if needed for Octave 3.4
    expected = expected';
end
t_ok(strcmp(got, expected), [t 'standard DC OPF']);
delete(tmpfnamedc);

if have_fcn('octave')
    warning(s1.state, 'Octave:load-file-in-path');
end
warning(s2.state, 'MATLAB:singularMatrix');

t_end;
