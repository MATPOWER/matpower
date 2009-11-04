function t_begin(num_of_tests, quiet)
%T_BEGIN  Begin running tests.
%   t_begin(num_of_tests, quiet) initializes the global test counters,
%   setting everything up to execute num_of_tests tests using t_ok()
%   and t_is(). If quiet is true, it will not print anything for the
%   individual tests, only a summary when t_end() is called.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

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
