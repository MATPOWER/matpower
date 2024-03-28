#####################
MP-Test User's Manual
#####################

..
    .. note::
    
       The new web-based version of the MP-Test User's Manual is not yet available. Please, continue to use the |MPTESTman| on GitHub for now.

|MPTEST>| is a set of functions for implementing unit testing in |MATLAB>| or
|Octave>|. It was initially developed for |MATPOWER>|, and is used by
|MATPOWER>|, |MIPSname| (|MIPS>|), |MPOM>|, and |MOSTname| (|MOST>|),
among many other packages. It also includes a function
:func:`have_feature` for detecting support for optional functionality.

Installation
============

**Note to** |*MATPOWER*| **users:** *MP-Test is included when you install*
|/MATPOWER/|. *There is generally no need to install it separately. You
can skip directly to step 3 to verify.*

Installation and use of |MPTEST>| requires familiarity with the basic operation
of MATLAB or Octave, including setting up your MATLAB/Octave path.

1.  Clone the repository or download and extract the zip file of the MP-Test
    distribution from the |MPTEST>| `project page
    <https://github.com/MATPOWER/mptest>`_ to the location of your
    choice. The files in the resulting ``mptest`` or ``mptestXXX`` directory,
    where ``XXX`` depends on the version of MP-Test, should not need to be
    modified, so it is recommended that they be kept separate from your
    own code. We will use :samp:`{<MPTEST>}` to denote the path to this
    directory.

2.  Add the following directories to your MATLAB/Octave path:

    - :samp:`{<MPTEST>}/lib`
    - :samp:`{<MPTEST>}/lib/t`

3.  At the |MATLAB|/Octave prompt, type ``test_mptest`` to run the test suite
    and verify that MP-Test is properly installed and functioning. The result
    should resemble the following::

        >> test_mptest
        t_test_fcns.......ok
        t_have_feature....ok
        All tests successful (29 of 29)
        Elapsed time 0.18 seconds.


Usage
=====

Write test functions of the following form, where :func:`t_ok`,
:func:`t_is`, :func:`t_str_match`, and :func:`t_file_match` are used to
test for specific conditions or matches.
::

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

Then run your test function:

.. code-block:: text

    >> mptest_ex1
    1..4
    ok 1 - size of pi
    skipped 2..3 - foo and bar tests require Octave
    ok 4 - 2+2 still equals 4
    All tests successful (2 passed, 2 skipped of 4)
    Elapsed time 0.00 seconds.

If you have several test functions, create a function to run them all as
follows::

    function test_everything_ex1(verbose)
    if nargin < 1
        verbose = 0;
    end
    tests = {};
    tests{end+1} = 'mptest_ex1';
    tests{end+1} = 't_test_fcns';
    
    t_run_tests( tests, verbose );

Run all of your tests at once. The output may look something like:

.. code-block:: text

    >> test_everything_ex1
    mptest_ex1.....ok (2 of 4 skipped)
    t_test_fcns....ok
    All tests successful (7 passed, 2 skipped of 9)
    Elapsed time 0.09 seconds.


Documentation
=============

The primary sources of documentation for MP-Test are this User's Manual, especially the :ref:`sec_mptest_reference` section, and the built-in ``help`` command. As with the built-in functions and toolbox routines in |MATLAB>| and |Octave>|, you can type ``help`` followed by the name of a command or M-file to get help on that particular function.

.. toctree::

   reference


Contributing
============

Please see our `contributing guidelines <https://github.com/MATPOWER/mptest/blob/master/CONTRIBUTING.md>`_ for details on how to contribute to the project or report issues.

License
=======

MP-Test is distributed under the `3-clause BSD license <https://github.com/MATPOWER/mptest/blob/master/LICENSE>`_.

Acknowledgments
===============

This material is based upon work supported in part by the |CERTS| and the Office of Electricity Delivery and Energy Reliability, Transmission Reliability Program of the U.S. Department of Energy under the National Energy Technology Laboratory Cooperative Agreement No. DE-FC26-09NT43321 and by the National Science Foundation under Grant Nos. 0532744, 1642341 and 1931421. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the funding agencies.
