function mptest_ex1(quiet)
if nargin < 1
quiet = 0;
end
t_begin(4, quiet);
t_ok(pi > 3, 'size of pi');
if exist('my_unimplemented_functionality', 'file')
t_ok(1, 'unimplemented_test1');
t_ok(1, 'unimplemented_test2');
else
t_skip(2, 'not yet written');
end
t_is(2+2, 4, 12, '2+2 still equals 4');
t_end;
