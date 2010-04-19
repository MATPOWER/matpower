function t_begin(num_of_tests, quiet)
%T_BEGIN  Begin running tests.
%   T_BEGIN(NUM_OF_TESTS, QUIET) initializes the global test counters,
%   setting everything up to execute NUM_OF_TESTS tests using T_OK
%   and T_IS. If QUIET is true, it will not print anything for the
%   individual tests, only a summary when T_END is called.
%
%   Example:
%       quiet = 0;
%       t_begin(5, quiet);
%       t_ok(pi > 3, 'size of pi');
%       t_skip(3, 'not yet written');
%       t_is(2+2, 4, 12, '2+2 still equals 4');
%       t_end;
%
%   See also T_END, T_OK, T_IS, T_SKIP, T_RUN_TESTS.

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

global t_quiet;
global t_num_of_tests;
global t_counter;
global t_ok_cnt;
global t_not_ok_cnt;
global t_skip_cnt;
global t_clock;

if nargin < 2
    quiet = 0;
end

t_quiet = quiet;
t_num_of_tests = num_of_tests;
t_counter = 1;
t_ok_cnt = 0;
t_not_ok_cnt = 0;
t_skip_cnt = 0;
t_clock = clock;

if ~t_quiet
    fprintf('1..%d\n', num_of_tests);
end
