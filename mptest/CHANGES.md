Change history for MP-Test
==========================


Version 7.0b1 - *Oct 30, 2016*
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
