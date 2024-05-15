Change history for MP-Test
==========================


Version 8.0 - *May 10, 2024*
----------------------------

#### 5/10/24
  - Release 8.0.

#### 3/1/24
  - Add Sphinx-based User and Reference documentation.


Version 8.0b1 - *Dec 12, 2022*
------------------------------

#### 12/12/22
  - Release 8.0b1.

#### 11/10/22
  - Add `t_str_match()` to test that a string/char array matches an
    expected value. Includes the option to apply a set of regular expression
    or simple string replacements before comparing.
  - Add `t_file_match()` to test that the contents of two text files
    match, with the option to delete one if they do match. Includes
    the option to apply a set of string or regular expresssion replacements
    before comparing.

#### 4/10/21
  - Allow `logical` type (i.e. `true` and `false`) as second input to
    `have_feature()`, as well as `numeric` type (1 and 0), to facilitate
    using (`logical`) output of prior call as toggle input.


Version 7.1 - *Oct 8, 2020*
---------------------------

#### 10/8/20
  - Release 7.1.

#### 9/22/20
  - Add `have_feature()`, moved from [MP-Opt-Model][2], as a modular,
    extensible alternative to `have_fcn()`, originally from [MATPOWER][3],
    where the detection of a feature named `<tag>` is implemented by the
    function `have_feature_<tag>()`. Includes feature detection functions
    for `'matlab'` and `'octave'`.
  - Switch from Travis-CI to GitHub workflows for continuous integration
    testing.

#### 9/21/20
  - Add `mptestver()` defining explicit version number.

#### 12/6/19
  - Improve handling of complex quantities in `t_is()`. Now displays
    differences in imaginary parts as well when there is a mismatch.
    Previously, it would show only the real parts which could be identical.


Version 7.0 - *Jun 20, 2019*
----------------------------

#### 6/20/19
  - Release 7.0.


Version 7.0b1 - *Oct 30, 2018*
------------------------------

#### 10/30/18
  - Release 7.0b1.
  - Fix bug in `t_is()` where comparing an integer value with a double
    value would pass when it should not.

#### 3/7/18
  - Replace `clock()`/`etime()` with `tic()`/`toc()` for timing.

#### 1/30/18
  - Added `abs()` to output of failed tests for consistency when
    comparing complex numbers (affects display only, test were correct).

#### 8/1/17
  - Modified `t_is()` to handle sparse matrix inputs.

#### 12/21/16
  - Updates for [Travis-CI][1] integration, with option to exit Octave
    or MATLAB if any test fails.
    *Thanks to Richard Lincoln for getting us started with Travis-CI.*


Version 6.0 - *Dec 16, 2016*
----------------------------

  - _no change_

#### 12/15/16
  - Moved development to GitHub: <https://github.com/MATPOWER/mptest>.



Version 6.0b2 - *Nov 1, 2016*
-----------------------------

  - _no change_


Version 6.0b2 - *Jun 1, 2016*
-----------------------------

#### 10/15/15
  - Modified `t_is()` to handle matrix inputs of dimension
    greater than two.
  - Added `t_test_fcns()` to test `t_ok()` and `t_is()` and manually
    check output of failed tests.


Version 5.1 - *Mar 20, 2015*
----------------------------

#### 2/25/15
  - Switch to more permissive 3-clause BSD license from GPL 3.0.


Version 5.0 - *Dec 17, 2014*
----------------------------

#### 8/11/14
  - Optionally return success flag from `t_ok()` and `t_is()`.


Version 5.0b1 - *Jul 1, 2014*
-----------------------------

#### 3/13/13
  - Empty `got` and `expected` arguments to `t_is()` now
    count as a passing test instead of an error, as long as
    the dimensions match.


Version 4.1 - *Dec 14, 2011*
-----------------------------

#### 6/29/11
  - Updated `t_is()` to properly print when result includes NaNs.


Version 4.0 - *Feb 7, 2011*
---------------------------

  - _no change_


Version 4.0b5 - *Dec 13, 2010*
------------------------------

#### 12/2/10
  - Improved output of `t_is()`. Includes only elements violating
    tolerance.


Version 4.0b4 - *May 21, 2010*
------------------------------

  - _no change_


Version 4.0b3 - *Apr 19, 2010*
------------------------------

#### 4/19/10
  - Changed licensing to GNU General Public license. See
    `LICENSE` and `COPYING` files for details.


Version 4.0b2 - *Mar 19, 2010*
------------------------------

#### 3/10/10
  - Massive help text update to more closely match MathWorks
    conventions; function names in ALL CAPS, See also ..., Examples, etc.


Version 4.0b1 - *Dec 24, 2009*
------------------------------

#### 11/4/09
  - Removed unnecessary `return` statement at end of all
    M-files. If anything it should be an `end` statement, but even
    that is optional, so we just let functions get terminated by the
    end-of-file or another function declaration.


Version 3.2 - *Sep 21, 2007*
----------------------------

  - _no change_


Version 3.1b2 - *Sep 15, 2006*
------------------------------

  - _no change_


Version 3.1b1 - *Aug 1, 2006*
-----------------------------

#### 10/14/05
  - Added total tests to printing of number of tests skipped.

#### 7/8/05
  - Updated `t_is()` to handle dimension mismatches.


Version 3.0 - *Feb 14, 2005*
----------------------------

  - _no change_


Version 3.0b4 - *Jan 28, 2005*
------------------------------

  - _no change_


Version 3.0b3 - *Sep 20, 2004*
------------------------------

  - _no change_


Version 3.0b2 - *Sep 7, 2004*
-----------------------------

#### 9/7/04
  - Modified `t_is()` to print max diff and tolerance upon failure.


Version 3.0b1 - *Aug 25, 2004*
------------------------------

#### 8/25/04
  - Added skipping of tests and reporting of skipped tests.

#### 7/15/04
  - Added `t` subdirectory with various tests and testing tools.

----
[1]: https://travis-ci.org
[2]: https://github.com/MATPOWER/mp-opt-model
[3]: https://github.com/MATPOWER/matpower
