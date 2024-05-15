What's New in MP-Test 8.0
-------------------------

#### Released May 10, 2024

Below is a summary of the changes since version 7.1 of MP-Test. See the
[`CHANGES.md`][1] file for all the gory details.

#### New Features:
  - New function `t_str_match()` to test that a string/char array matches an
    expected value. Includes the option to apply a set of regular expression
    or simple string replacements before comparing.
  - New function `t_file_match()` to test that the contents of two text files
    match, with the option to delete one if they do match. Includes
    the option to apply a set of string or regular expresssion replacements
    before comparing.

#### New Documentation:
  - New Sphinx-based [User's Manual][2] and [Reference documentation][3].

#### Other Change:
  - Allow `logical` type (i.e. `true` and `false`) as second input to
    `have_feature()`, as well as `numeric` type (1 and 0), to facilitate
    using (`logical`) output of prior call as toggle input.


[1]: ../../CHANGES.md
[2]: https://matpower.org/doc/mptest/
[3]: https://matpower.org/doc/mptest/reference.html
