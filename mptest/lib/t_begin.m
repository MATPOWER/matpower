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

%   MP-Test
%   Copyright (c) 2004-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

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
t_clock = tic;

if ~t_quiet
    fprintf('1..%d\n', num_of_tests);
end
