function test_everything_ex1(verbose)
if nargin < 1
verbose = 0;
end
tests = {};
tests{end+1} = 'mptest_ex1';
tests{end+1} = 't_test_fcns';

t_run_tests( tests, verbose );
