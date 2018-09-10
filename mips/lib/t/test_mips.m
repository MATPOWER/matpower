function success = test_mips(verbose, exit_on_fail)
%TEST_MIPS  Run all MIPS tests.
%   TEST_MIPS
%   TEST_MIPS(VERBOSE)
%   TEST_MIPS(VERBOSE, EXIT_ON_FAIL)
%   SUCCESS = TEST_MIPS(...)
%
%   Runs all of the MIPS tests. If VERBOSE is true (false by default),
%   it prints the details of the individual tests. If EXIT_ON_FAIL is true
%   (false by default), it will exit MATLAB or Octave with a status of 1
%   unless T_RUN_TESTS returns ALL_OK.
%
%   See also T_RUN_TESTS.

%   MIPS
%   Copyright (c) 2016-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MIPS.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mips for more info.

if nargin < 2
    exit_on_fail = 0;
    if nargin < 1
        verbose = 0;
    end
end

tests = {};

%% MIPS tests
tests{end+1} = 't_mplinsolve';
tests{end+1} = 't_mips';
tests{end+1} = 't_mips_pardiso';
tests{end+1} = 't_qps_mips';

%% run the tests
all_ok = t_run_tests( tests, verbose );

%% handle success/failure
if exit_on_fail && ~all_ok
    exit(1);
end
if nargout
    success = all_ok;
end
