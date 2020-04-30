function success = test_mp_opt_model(verbose, exit_on_fail)
%TEST_MP_OPT_MODEL  Run all MATPOWER tests.
%   TEST_MP_OPT_MODEL
%   TEST_MP_OPT_MODEL(VERBOSE)
%   TEST_MP_OPT_MODEL(VERBOSE, EXIT_ON_FAIL)
%   SUCCESS = TEST_MP_OPT_MODEL(...)
%
%   Runs all of the MP-Opt-Model tests. If VERBOSE is true (false by default),
%   it prints the details of the individual tests. If EXIT_ON_FAIL is true
%   (false by default), it will exit MATLAB or Octave with a status of 1
%   unless T_RUN_TESTS returns ALL_OK.
%
%   See also T_RUN_TESTS.

%   MATPOWER
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 2
    exit_on_fail = 0;
    if nargin < 1
        verbose = 0;
    end
end

tests = {};

%% MATPOWER base test
tests{end+1} = 't_nested_struct_copy';
tests{end+1} = 't_have_fcn';
tests{end+1} = 't_mips';
tests{end+1} = 't_mips_pardiso';
tests{end+1} = 't_qps_mips';
tests{end+1} = 't_qps_matpower';
tests{end+1} = 't_miqps_matpower';
tests{end+1} = 't_nlps_matpower';
tests{end+1} = 't_opt_model';

%% run the tests
all_ok = t_run_tests( tests, verbose );

%% handle success/failure
if exit_on_fail && ~all_ok
    exit(1);
end
if nargout
    success = all_ok;
end
