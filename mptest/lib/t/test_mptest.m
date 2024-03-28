function success = test_mptest(verbose, exit_on_fail)
%test_mptest - Run all MP-Test tests.
% ::
%
%   test_mptest
%   test_mptest(verbose)
%   test_mptest(verbose, exit_on_fail)
%   success = test_mptest(...)
%
% Runs all of the MP-Test tests. If ``verbose`` is true *(false by default)*,
% it prints the details of the individual tests. If ``exit_on_fail`` is true
% *(false by default)*, it will exit MATLAB or Octave with a status of 1
% unless t_run_tests returns ``all_ok`` true.
%
% See also t_run_tests.

%   MP-Test
%   Copyright (c) 2016-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

if nargin < 2
    exit_on_fail = 0;
    if nargin < 1
        verbose = 0;
    end
end

tests = {};

%% MP-Test base test
tests{end+1} = 't_test_fcns';
tests{end+1} = 't_have_feature';

%% run the tests
all_ok = t_run_tests( tests, verbose );

%% handle success/failure
if exit_on_fail && ~all_ok
    exit(1);
end
if nargout
    success = all_ok;
end
