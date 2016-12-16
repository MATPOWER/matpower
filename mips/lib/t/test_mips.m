function test_mips(verbose)
%TEST_MIPS  Run all MIPS tests.
%   TEST_MIPS runs all of the MIPS tests.
%   TEST_MIPS(VERBOSE) prints the details of the individual tests
%   if VERBOSE is true.
%
%   See also T_RUN_TESTS.

%   MIPS
%   Copyright (c) 2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MIPS.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mips for more info.

if nargin < 1
    verbose = 0;
end

tests = {};

%% MATPOWER base test
tests{end+1} = 't_mplinsolve';
tests{end+1} = 't_mips';
tests{end+1} = 't_qps_mips';

t_run_tests( tests, verbose );
