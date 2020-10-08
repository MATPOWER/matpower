What's New in MP-Test 7.1
-------------------------

#### Released Oct 8, 2020

Below is a summary of the changes since version 7.0 of MP-Test. See the
[`CHANGES.md`][1] file for all the gory details.

#### New Features:
  - New function `have_feature()` to detect availability and version
    information for optional functionality. This is a modular, extensible
    replacement for `have_fcn()`, from [MATPOWER][2] and [MP-Opt-Model][3].
    The detection of a feature named `<tag>` is implemented by the
    function `have_feature_<tag>()`. Includes feature detection functions
    for `'matlab'` and `'octave'`.

#### Other Changes:
  - Expanded documentation in [README.md][4].
  - Improve handling of complex quantities in `t_is()`. Now displays
    differences in imaginary parts as well when there is a mismatch.
    Previously, it would show only the real parts which could be identical.
  - Add `mptestver()` defining explicit version number.
  - Switch from Travis-CI to GitHub workflows for continuous integration
    testing.


[1]: ../../CHANGES.md
[2]: https://github.com/MATPOWER/matpower
[3]: https://github.com/MATPOWER/mp-opt-model
[4]: ../../README.md
