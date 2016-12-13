function test_most(verbose)
%TEST_MOST  Run all MOST tests.
%   TEST_MOST runs all of the MOST tests.
%   TEST_MOST(VERBOSE) prints the details of the individual tests
%   if VERBOSE is true.
%
%   See also T_RUN_TESTS.

%   MOST
%   Copyright (c) 2004-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    verbose = 0;
end

tests = {};

%% MATPOWER base test
have_c3sopf = exist('c3sopf', 'file') == 2;
tests{end+1} = 't_most_3b_1_1_0';
tests{end+1} = 't_most_3b_3_1_0';
if have_c3sopf
    tests{end+1} = 't_most_3b_1_1_2';
    tests{end+1} = 't_most_3b_3_1_2';
end
tests{end+1} = 't_most_30b_1_1_0';
tests{end+1} = 't_most_30b_3_1_0';
if have_c3sopf
    tests{end+1} = 't_most_30b_1_1_17';
    tests{end+1} = 't_most_30b_3_1_17';
end
tests{end+1} = 't_most_fixed_res';
tests{end+1} = 't_most_w_ds';
tests{end+1} = 't_most_30b_1_1_0_uc';
if have_c3sopf
    tests{end+1} = 't_most_sp';
    tests{end+1} = 't_most_spuc';
end
tests{end+1} = 't_most_uc';
tests{end+1} = 't_most_suc';

t_run_tests( tests, verbose );
