function t_begin(num_of_tests)

global t_num_of_tests;
global t_counter;
global t_ok_cnt;
global t_not_ok_cnt;
global t_clock;

fprintf('1..%d\n', num_of_tests);

t_num_of_tests = num_of_tests;
t_counter = 1;
t_ok_cnt = 0;
t_not_ok_cnt = 0;
t_clock = clock;

return;
