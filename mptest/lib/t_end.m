function t_begin(num_of_tests)

global t_num_of_tests;
global t_counter;
global t_ok_cnt;
global t_not_ok_cnt;

t_counter = t_counter - 1;

if t_counter == t_num_of_tests & t_ok_cnt == t_counter & t_not_ok_cnt == 0
	fprintf('All tests successful (%d of %d).\n', t_ok_cnt, t_num_of_tests);
else
	fprintf('Ran %d of %d tests: %d passed, %d failed\n', t_counter, t_num_of_tests, t_ok_cnt, t_not_ok_cnt);
end

return;
