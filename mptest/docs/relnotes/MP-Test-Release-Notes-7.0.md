What's New in MP-Test 7.0
-------------------------

#### Released Jun 20, 2019

Below is a summary of the changes since version 6.0 of MP-Test. See the
[`CHANGES.md`][1] file for all the gory details.

#### New Features:
  - Updates for [Travis-CI][2] integration, with option to exit Octave
    or MATLAB if any test fails.
    *Thanks to Richard Lincoln for getting us started with Travis-CI.*
  - Modified `t_is()` to handle sparse matrix inputs.

#### Bugs Fixed:
  - Fix bug in `t_is()` where comparing an integer value with a double
    value would pass when it should not.
  - Added `abs()` to output of failed tests for consistency when
    comparing complex numbers (affects display only, test were correct).

#### Other Changes:
  - Replace `clock()`/`etime()` with `tic()`/`toc()` for timing.


[1]: https://github.com/MATPOWER/mptest/blob/master/CHANGES.md
[2]: https://travis-ci.org
