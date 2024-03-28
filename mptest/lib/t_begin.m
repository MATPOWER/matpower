function t_begin(num_of_tests, quiet)
% t_begin - Begin running tests.
% ::
%
%   t_begin(num_of_tests)
%   t_begin(num_of_tests, quiet)
%
% Initializes the global test counters, setting everything up to execute
% ``num_of_tests`` tests using the individual test functions t_ok, t_is,
% t_file_match, t_str_match and t_skip. If ``quiet`` is true, it will not
% print anything for the individual tests, only a summary when t_end is
% called.
%
% Inputs:
%   num_of_tests (integer) : number of tests expected
%   quiet (boolean) : *(optional, default = false)* if true, prevents
%       printing individual tests results
%
% Example::
%
%   quiet = 0;
%   t_begin(5, quiet);
%   t_ok(pi > 3, 'size of pi');
%   t_skip(3, 'not yet written');
%   t_is(2+2, 4, 12, '2+2 still equals 4');
%   t_end;
%
% See also t_end, t_ok, t_is, t_str_match, t_file_match, t_skip, t_run_tests.

%   MP-Test
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
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
