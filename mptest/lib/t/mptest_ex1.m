function mptest_ex1(quiet)
if nargin < 1
    quiet = 0;
end
t_begin(4, quiet);
t_ok(pi > 3, 'size of pi');
if have_feature('octave')
    t_ok(1, 'Octave-only test foo');
    t_ok(1, 'Octave-only test bar');
else
    t_skip(2, 'foo and bar tests require Octave');
end
t_is(2+2, 4, 12, '2+2 still equals 4');
t_end;
