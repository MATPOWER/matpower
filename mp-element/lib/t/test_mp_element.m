function success = test_mp_element(verbose, exit_on_fail)
%T_MP_ELEMENT  Run all MATPOWER tests.
%   T_MP_ELEMENT
%   T_MP_ELEMENT(VERBOSE)
%   T_MP_ELEMENT(VERBOSE, EXIT_ON_FAIL)
%   SUCCESS = T_MP_ELEMENT(...)
%
%   Runs all of the MP-Element tests. If VERBOSE is true (false by default),
%   it prints the details of the individual tests. If EXIT_ON_FAIL is true
%   (false by default), it will exit MATLAB or Octave with a status of 1
%   unless T_RUN_TESTS returns ALL_OK.
%
%   See also T_RUN_TESTS.

%   MATPOWER
%   Copyright (c) 2004-2019, Power Systems Engineering Research Center (PSERC)
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
tests{end+1} = 't_mp_mapped_array';
tests{end+1} = 't_mp_table';
tests{end+1} = 't_mp_data_model';
tests{end+1} = 't_dmc_element';
tests{end+1} = 't_mp_dm_converter_mpc2';
tests{end+1} = 't_nm_element';
tests{end+1} = 't_port_inj_current_acc';
tests{end+1} = 't_port_inj_current_acp';
tests{end+1} = 't_port_inj_power_acc';
tests{end+1} = 't_port_inj_power_acp';
tests{end+1} = 't_node_test';
tests{end+1} = 't_run_mp';
tests{end+1} = 't_run_mp_3p';
tests{end+1} = 't_pretty_print';
tests{end+1} = 't_mpxt_reserves';

%% run the tests
all_ok = t_run_tests( tests, verbose );

%% handle success/failure
if exit_on_fail && ~all_ok
    exit(1);
end
if nargout
    success = all_ok;
end
