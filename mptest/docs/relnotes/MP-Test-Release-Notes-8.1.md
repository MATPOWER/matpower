What's New in MP-Test 8.1
-------------------------

#### Released July 5, 2025

Below is a summary of the changes since version 8.0 of MP-Test. See the
[`CHANGES.md`][1] file for all the gory details.

#### New Features:
  - Two new functions to assist with code for debugging:
    - `assert_debug()` -- calls `assert()` if `DEBUG_MODE` is on
    - `toggle_debug_mode()` -- set/toggles whether `DEBUG_MODE` is on

#### Changes:
  - Enhance `t_is()` to handle `Inf` values appropriately.


[1]: ../../CHANGES.md
[2]: https://matpower.org/doc/mptest/
[3]: https://matpower.org/doc/mptest/reference.html
