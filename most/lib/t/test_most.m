function test_mops(verbose)
%TEST_MOPS  Run all MOPS tests.
%   TEST_MOPS runs all of the MOPS tests.
%   TEST_MOPS(VERBOSE) prints the details of the individual tests
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
tests{end+1} = 't_apply_changes';
tests{end+1} = 't_mops_3b_1_1_0';
tests{end+1} = 't_mops_3b_3_1_0';
tests{end+1} = 't_mops_3b_1_1_2';
tests{end+1} = 't_mops_3b_3_1_2';
tests{end+1} = 't_mops_30b_1_1_0';
tests{end+1} = 't_mops_30b_3_1_0';
tests{end+1} = 't_mops_30b_1_1_17';
tests{end+1} = 't_mops_30b_3_1_17';
tests{end+1} = 't_mops_fixed_res';
tests{end+1} = 't_mops_w_ds';
tests{end+1} = 't_mops_30b_1_1_0_uc';
tests{end+1} = 't_mops_sp';
tests{end+1} = 't_mops_spuc';
tests{end+1} = 't_mops_uc';
tests{end+1} = 't_mops_suc';

t_run_tests( tests, verbose );
