MP-Test
=======

MP-Test is a set of functions for implementing unit testing in MATLAB or
Octave. It was initially developed for [MATPOWER][1], and is used by
[MATPOWER][1], [MATPOWER Interior Point Solver (MIPS)][2] and
[MATPOWER Optimal Scheduling Tool (MOST)][3].

Installation
------------

Installation and use of MP-Test requires familiarity with the basic operation
of MATLAB or Octave, including setting up your MATLAB path.

1.  Clone the repository or download and extract the zip file of the MIPS
    distribution from the [MP-Test project page][4] to the location of your
    choice. The files in the resulting `mptest` or `mptestXXX` directory,
    where `XXX` depends on the version of MP-Test, should not need to be
    modified, so it is recommended that they be kept separate from your
    own code. We will use `<MPTEST>` to denote the path to this directory.

2.  Add the following directories to your MATLAB/Octave path:
    *   `<MPTEST>/lib`
    *   `<MPTEST>/lib/t`

3.  At the MATLAB prompt, type `test_mptest` to run the test suite and
    verify that MP-Test is properly installed and functioning. The result
    should resemble the following:
```matlab
  >> test_mptest
  t_test_fcns....ok
  All tests successful (1 of 1)
  Elapsed time 0.01 seconds.
```

Usage
-----

*   Write test functions of the following form, where `t_ok` and `t_is` are
    used to test for specific conditions or matches, respectively.
```matlab
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
```

*   Then run your test function:
```
  >> mptest_ex1
  1..4
  ok 1 - size of pi
  skipped tests 2..3 : not yet written
  ok 4 - 2+2 still equals 4
  All tests successful (2 passed, 2 skipped of 4)
  Elapsed time 0.00 seconds.
```

*   If you have several test functions, create a function to run them all as follows:
```matlab
  function test_everything_ex1(verbose)
  if nargin < 1
    verbose = 0;
  end
  tests = {};
  tests{end+1} = 'mptest_ex1';
  tests{end+1} = 't_test_fcns';

  t_run_tests( tests, verbose );
```

*   Run all of your tests at once. The output may look something like:
```
  >> test_everything_ex1
  mptest_ex1.....ok (2 of 4 skipped)
  t_test_fcns....ok
  All tests successful (3 passed, 2 skipped of 5)
  Elapsed time 0.00 seconds.
```

Documentation
-------------

The primary source of documentation for MP-Test is the built-in `help`
command. As with the built-in functions and toolbox routines in MATLAB
and Octave, you can type `help` followed by the name of a command or
M-file to get help on that particular function. All of the M-files in
MP-Test have such documentation and this should be considered the main
reference for the calling options for each function, namely:: `t_begin`, `t_end`, `t_ok`, `t_is`, `t_skip` and `t_run_tests`.

Contributing
------------

Please see our [contributing guidelines][5] for details on how to
contribute to the project or report issues.

License
-------

MP-Test is distributed under the [3-clause BSD license][6].

----
[1]: https://github.com/MATPOWER/matpower
[2]: https://github.com/MATPOWER/mips
[3]: https://github.com/MATPOWER/most
[4]: https://github.com/MATPOWER/mptest
[5]: CONTRIBUTING.md
[6]: LICENSE