function success = test_most(verbose, exit_on_fail)
% test_most - Run all MOST tests.
% ::
%
%   test_most
%   test_most(verbose)
%   test_most(verbose, exit_on_fail)
%   success = test_most(...)
%
% Runs all of the MOST tests. If ``verbose`` is true *(false by default)*,
% it prints the details of the individual tests. If ``exit_on_fail`` is true
% *(false by default)*, it will exit MATLAB or Octave with a status of 1
% unless t_run_tests returns ``all_ok`` true.
%
% See also t_run_tests.

%   MOST
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

if nargin < 2
    exit_on_fail = 0;
    if nargin < 1
        verbose = 0;
    end
end

tests = {};

%% MOST tests
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
tests{end+1} = 't_most_30b_1_1_0_uc';
if have_c3sopf
    tests{end+1} = 't_most_sp';
    tests{end+1} = 't_most_spuc';
end
tests{end+1} = 't_most_mpopf';
tests{end+1} = 't_most_uc';
tests{end+1} = 't_most_suc';
tests{end+1} = 't_most_tlmp';
tests{end+1} = 't_most_w_ds';

%% run the tests
all_ok = t_run_tests( tests, verbose );

%% handle success/failure
if exit_on_fail && ~all_ok
    exit(1);
end
if nargout
    success = all_ok;
end
