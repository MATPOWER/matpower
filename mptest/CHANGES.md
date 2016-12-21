Change history for MP-Test
==========================


Since last release
------------------

* 2016-12-21 - Updates for [Travis-CI][1] integration (thanks to Richard
  Lincoln), with option to exit Octave or Matlab if any test fails.


6.0 (2016-12-16)
----------------

* _no change_


6.0b2 (2016-11-01)
------------------

* _no change_


6.0b1 (2016-06-01)
------------------

* 2015-10-15 - Modified `t_is()` to handle matrix inputs of dimension
  greater than two.

* 2015-10-15 - Added `t_test_fcns()` to test `t_ok()` and `t_is()` and manually
  check output of failed tests.


5.1 (2015-03-20)
----------------

* 2015-02-25 - Switch to more permissive 3-clause BSD license from GPL 3.0.


5.0 (2014-12-17)
----------------

* 2014-08-11 - Optionally return success flag from `t_ok()` and `t_is()`.


5.0b1 (2014-07-01)
------------------

* 2013-03-13 - Empty `got` and `expected` arguments to `t_is()` now
  count as a passing test instead of an error, as long as
  the dimensions match.


4.1 (2011-12-14)
----------------

* 2011-06-29 - Updated `t_is()` to properly print when result includes NaNs.


4.0 (2011-02-07)
----------------

* _no change_


4.0b5 (2010-12-13)
------------------

* 2010-12-02 - Improved output of `t_is()`. Includes only elements violating
  tolerance.


4.0b4 (2010-05-21)
------------------

* _no change_


4.0b3 (2010-04-19)
------------------

* 2010-04-19 - Changed licensing to GNU General Public license. See
  `LICENSE` and `COPYING` files for details.


4.0b2 (2010-03-19)
------------------

* 2010-03-10 - Massive help text update to more closely match MathWorks
  conventions; function names in ALL CAPS, See also ..., Examples, etc.


4.0b1 (2009-12-24)
-----------------

* 2009-11-04 - Removed unnecessary `return` statement at end of all
  M-files. If anything it should be an `end` statement, but even
  that is optional, so we just let functions get terminated by the
  end-of-file or another function declaration.


3.2 (2007-09-21)
----------------

* _no change_


3.1b2 (2006-09-15)
------------------

* _no change_


3.1b1 (2006-08-01)
------------------

* 2005-10-14 - Added total tests to printing of number of tests skipped.

* 2005-07-08 - Updated `t_is()` to handle dimension mismatches.
    

3.0 (2005-02-14)
----------------

* _no change_


3.0b4 (2005-01-28)
------------------

* _no change_


3.0b3 (2004-09-20)
-------------------

* _no change_


3.0b2 (2004-09-07)
-------------------

* 2004-09-07 - Modified `t_is()` to print max diff and tolerance upon failure.


3.0b1 (2004-08-25)
-----------------

* 2004-08-25 - Added skipping of tests and reporting of skipped tests.

* 2004-07-15 - Added `t` subdirectory with various tests and testing tools.

----
[1]: https://travis-ci.org
