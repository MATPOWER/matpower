function all_ok_ = t_run_tests(test_names, verbose)
% t_run_tests - Run a series of tests.
% ::
%
%   all_ok = t_run_tests(test_names, verbose)
%
% Inputs:
%   test_names (cell array) : cell array with names of individual test
%       scripts to run
%   verbose (integer) : *(optional, default = 0)* level of detail of output
%
%       - 0 -- print a single line for each test, plus summary line
%       - 1 -- print full output for each test, plus summary section
%
% Output:
%   all_ok (boolean) : *(optional)* true if all tests passed and the number
%       of tests matches the expected number, false otherwise
%
% Runs a set of tests whose names are given in the cell array ``test_names``.
%
% Example::
%
%   tests{end+1} = 't_loadcase';
%   tests{end+1} = 't_jacobian';
%   tests{end+1} = 't_hessian';
%   t_run_tests( tests, verbose );
%
% See also t_begin, t_end.

%   MP-Test
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

if nargin < 2
    verbose = 0;
end

global t_num_of_tests;
global t_counter;
global t_ok_cnt;
global t_not_ok_cnt;
global t_skip_cnt;

%% figure out padding for printing
if ~verbose
    len = zeros(length(test_names), 1);
    for k = 1:length(test_names)
        len(k) = length(test_names{k});
    end
    maxlen = max(len);
end

%% initialize statistics
num_of_tests = 0;
counter = 0;
ok_cnt = 0;
not_ok_cnt = 0;
skip_cnt = 0;

t0 = tic;
for k = 1:length(test_names)
    if verbose
        fprintf('\n----------  %s  ----------\n', test_names{k});
    else
        pad = maxlen + 4 - length(test_names{k});
        fprintf('%s', test_names{k});
        for m = 1:pad, fprintf('.'); end
    end
    feval( test_names{k}, ~verbose );
    
    num_of_tests    = num_of_tests  + t_num_of_tests;
    counter         = counter       + t_counter;
    ok_cnt          = ok_cnt        + t_ok_cnt;
    not_ok_cnt      = not_ok_cnt    + t_not_ok_cnt;
    skip_cnt        = skip_cnt      + t_skip_cnt;
end

if verbose
    fprintf('\n\n----------  Summary  ----------\n');
end
if counter == num_of_tests && counter == ok_cnt + skip_cnt && not_ok_cnt == 0
    all_ok = 1;
    if skip_cnt
        fprintf('All tests successful (%d passed, %d skipped of %d)', ...
            ok_cnt, skip_cnt, num_of_tests);
    else
        fprintf('All tests successful (%d of %d)', ok_cnt, num_of_tests);
    end
else
    all_ok = 0;
    fprintf('Ran %d of %d tests: %d passed, %d failed', ...
        counter, num_of_tests, ok_cnt, not_ok_cnt);
    if skip_cnt
        fprintf(', %d skipped', skip_cnt);
    end
end
fprintf('\nElapsed time %.2f seconds.\n', toc(t0));

if nargout
    all_ok_ = all_ok;   %% copy to optional output arg
end
