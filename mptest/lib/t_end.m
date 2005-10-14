function t_end
%T_END  Finish running tests and print statistics.
%   t_end checks the global counters that were updated by calls to
%   t_ok() and t_is() and prints out a summary of the test results.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

global t_quiet;
global t_num_of_tests;
global t_counter;
global t_ok_cnt;
global t_not_ok_cnt;
global t_skip_cnt;
global t_clock;

t_counter = t_counter - 1;

if t_counter == t_num_of_tests & ...
        t_counter == t_ok_cnt + t_skip_cnt & ...
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
    fprintf('\nElapsed time %.2f seconds.\n', etime(clock, t_clock));
end

return;
