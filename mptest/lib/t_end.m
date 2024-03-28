function t_end
% t_end - Finish running tests and print statistics.
% ::
%
%   t_end
%
% Checks the global counters that were updated by calls to the individual
% test functions t_ok, t_is, t_file_match, t_str_match and t_skip, and
% prints out a summary of the test results.
%
% Example::
%
%   quiet = 0;
%   t_begin(5, quiet);
%   t_ok(pi > 3, 'size of pi');
%   t_skip(3, 'not yet written');
%   t_is(2+2, 4, 12, '2+2 still equals 4');
%   t_end;
%
% See also t_begin, t_ok, t_is, t_str_match, t_file_match, t_skip, t_run_tests.

%   MP-Test
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

global t_quiet;
global t_num_of_tests;
global t_counter;
global t_ok_cnt;
global t_not_ok_cnt;
global t_skip_cnt;
global t_clock;

t_counter = t_counter - 1;

if t_counter == t_num_of_tests && ...
        t_counter == t_ok_cnt + t_skip_cnt && ...
        t_not_ok_cnt == 0
    all_ok = 1;
else
    all_ok = 0;
end

if t_quiet
    if all_ok
        fprintf('ok');
        if t_skip_cnt
            fprintf(' (%d of %d skipped)', t_skip_cnt, t_num_of_tests);
        end
    else
        fprintf('not ok\n');
        fprintf('\t#####  Ran %d of %d tests: %d passed, %d failed', ...
            t_counter, t_num_of_tests, t_ok_cnt, t_not_ok_cnt);
        if t_skip_cnt
            fprintf(', %d skipped', t_skip_cnt);
        end
    end
    fprintf('\n');
else
    if all_ok
        if t_skip_cnt
            fprintf('All tests successful (%d passed, %d skipped of %d)', ...
                t_ok_cnt, t_skip_cnt, t_num_of_tests);
        else
            fprintf('All tests successful (%d of %d)', t_ok_cnt, t_num_of_tests);
        end
    else
        fprintf('Ran %d of %d tests: %d passed, %d failed', ...
            t_counter, t_num_of_tests, t_ok_cnt, t_not_ok_cnt);
        if t_skip_cnt
            fprintf(', %d skipped', t_skip_cnt);
        end
    end
    fprintf('\nElapsed time %.2f seconds.\n', toc(t_clock));
end
