function test_most(verbose)
%TEST_MOST  Run all MOST tests.
%   TEST_MOST runs all of the MOST tests.
%   TEST_MOST(VERBOSE) prints the details of the individual tests
%   if VERBOSE is true.
%
%   See also T_RUN_TESTS.

%   MATPOWER
%   Copyright (c) 2004-2015 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   $Id: test_matpower.m 2703 2015-04-14 15:39:38Z ray $
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    verbose = 0;
end

tests = {};

%% MATPOWER base test
tests{end+1} = 't_most_3b_1_1_0';
tests{end+1} = 't_most_3b_3_1_0';
tests{end+1} = 't_most_3b_1_1_2';
tests{end+1} = 't_most_3b_3_1_2';
tests{end+1} = 't_most_30b_1_1_0';
tests{end+1} = 't_most_30b_3_1_0';
tests{end+1} = 't_most_30b_1_1_17';
tests{end+1} = 't_most_30b_3_1_17';
tests{end+1} = 't_most_fixed_res';
tests{end+1} = 't_most_w_ds';
tests{end+1} = 't_most_30b_1_1_0_uc';
tests{end+1} = 't_most_sp';
tests{end+1} = 't_most_spuc';
tests{end+1} = 't_most_uc';
tests{end+1} = 't_most_suc';

t_run_tests( tests, verbose );
