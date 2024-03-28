function success = test_mp_opt_model(verbose, exit_on_fail)
% test_mp_opt_model - Run all MP-Opt-Model tests.
% ::
%
%   test_mp_opt_model
%   test_mp_opt_model(verbose)
%   test_mp_opt_model(verbose, exit_on_fail)
%   success = test_mp_opt_model(...)
%
% Runs all of the MP-Opt-Model tests. If ``verbose`` is true *(false by
% default)*, it prints the details of the individual tests. If
% ``exit_on_fail`` is true *(false by default)*, it will exit MATLAB or
% Octave with a status of 1 unless t_run_tests returns ``all_ok`` true.
%
% See also t_run_tests.

%   MP-Opt-Model
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin < 2
    exit_on_fail = 0;
    if nargin < 1
        verbose = 0;
    end
end

tests = {};

%% MP-Opt-Model tests
tests{end+1} = 't_have_fcn';
tests{end+1} = 't_nested_struct_copy';
tests{end+1} = 't_nleqs_master';
tests{end+1} = 't_pnes_master';
tests{end+1} = 't_qps_master';
tests{end+1} = 't_miqps_master';
tests{end+1} = 't_nlps_master';
tests{end+1} = 't_opt_model';
tests{end+1} = 't_om_solve_leqs';
tests{end+1} = 't_om_solve_nleqs';
tests{end+1} = 't_om_solve_pne';
tests{end+1} = 't_om_solve_qps';
tests{end+1} = 't_om_solve_miqps';
tests{end+1} = 't_om_solve_nlps';

%% run the tests
all_ok = t_run_tests( tests, verbose );

%% handle success/failure
if exit_on_fail && ~all_ok
    exit(1);
end
if nargout
    success = all_ok;
end
