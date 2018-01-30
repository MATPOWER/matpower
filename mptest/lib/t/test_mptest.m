function success = test_mptest(verbose, exit_on_fail)
%TEST_MPTEST  Run all MPTEST tests.
%   TEST_MPTEST
%   TEST_MPTEST(VERBOSE)
%   TEST_MPTEST(VERBOSE, EXIT_ON_FAIL)
%   SUCCESS = TEST_MPTEST(...)
%
%   Runs all of the MPTEST tests. If VERBOSE is true (false by default),
%   it prints the details of the individual tests. If EXIT_ON_FAIL is true
%   (false by default), it will exit MATLAB or Octave with a status of 1
%   unless T_RUN_TESTS returns ALL_OK.
%
%   See also T_RUN_TESTS.

%   MP-Test
%   Copyright (c) 2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

if nargin < 2
    exit_on_fail = 0;
    if nargin < 1
        verbose = 0;
    end
end

tests = {};

%% MP-Test base test
tests{end+1} = 't_test_fcns';

%% run the tests
all_ok = t_run_tests( tests, verbose );

%% handle success/failure
if exit_on_fail && ~all_ok
    exit(1);
end
if nargout
    success = all_ok;
end
