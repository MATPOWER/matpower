function test_mptest(verbose)
%TEST_MPTEST  Run all MPTEST tests.
%   TEST_MPTEST runs all of the MPTEST tests.
%   TEST_MPTEST(VERBOSE) prints the details of the individual tests
%   if VERBOSE is true.
%
%   See also T_RUN_TESTS.

%   MP-Test
%   Copyright (c) 2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

if nargin < 1
    verbose = 0;
end

tests = {};

%% MP-Test base test
tests{end+1} = 't_test_fcns';

t_run_tests( tests, verbose );
