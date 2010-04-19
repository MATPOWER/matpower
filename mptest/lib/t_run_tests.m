function t_run_tests(test_names, verbose)
%T_RUN_TESTS  Run a series of tests.
%   T_RUN_TESTS(TEST_NAMES, VERBOSE) runs a set of tests whose names
%   are given in the cell array TEST_NAMES. If the optional parameter
%   VERBOSE is true, it prints the details of the individual tests.
%
%   Example:
%       tests{end+1} = 't_loadcase';
%       tests{end+1} = 't_jacobian';
%       tests{end+1} = 't_hessian';
%       t_run_tests( tests, verbose );
%
%   See also T_BEGIN, T_END.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2004-2010 by Power System Engineering Research Center (PSERC)
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
%   other modules (such as M-files and MEX-files) available in a
%   Matlab (or compatible) environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

if nargin < 2
    verbose = 0;
end

global t_num_of_tests;
global t_counter;
global t_ok_cnt;
global t_not_ok_cnt;
global t_skip_cnt;

%% figure out padding for printing
if ~verbose
    len = zeros(length(test_names), 1);
    for k = 1:length(test_names)
        len(k) = length(test_names{k});
    end
    maxlen = max(len);
end

%% initialize statistics
num_of_tests = 0;
counter = 0;
ok_cnt = 0;
not_ok_cnt = 0;
skip_cnt = 0;

t0 = clock;
for k = 1:length(test_names)
    if verbose
        fprintf('\n----------  %s  ----------\n', test_names{k});
    else
        pad = maxlen + 4 - length(test_names{k});
        fprintf('%s', test_names{k});
        for m = 1:pad, fprintf('.'); end
    end
    feval( test_names{k}, ~verbose );
    
    num_of_tests    = num_of_tests  + t_num_of_tests;
    counter         = counter       + t_counter;
    ok_cnt          = ok_cnt        + t_ok_cnt;
    not_ok_cnt      = not_ok_cnt    + t_not_ok_cnt;
    skip_cnt        = skip_cnt      + t_skip_cnt;
end

if verbose
    fprintf('\n\n----------  Summary  ----------\n');
end
if counter == num_of_tests && counter == ok_cnt + skip_cnt && not_ok_cnt == 0
    if skip_cnt
        fprintf('All tests successful (%d passed, %d skipped of %d)', ...
            ok_cnt, skip_cnt, num_of_tests);
    else
        fprintf('All tests successful (%d of %d)', ok_cnt, num_of_tests);
    end
else
    fprintf('Ran %d of %d tests: %d passed, %d failed', ...
        counter, num_of_tests, ok_cnt, not_ok_cnt);
    if skip_cnt
        fprintf(', %d skipped', skip_cnt);
    end
end
fprintf('\nElapsed time %.2f seconds.\n', etime(clock, t0));
