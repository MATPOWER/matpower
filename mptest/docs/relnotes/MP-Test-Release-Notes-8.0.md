What's New in MP-Test 8.0b1
---------------------------

#### Released Dec 12, 2022

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

#### Other Change:
  - Allow `logical` type (i.e. `true` and `false`) as second input to
    `have_feature()`, as well as `numeric` type (1 and 0), to facilitate
    using (`logical`) output of prior call as toggle input.


[1]: ../../CHANGES.md
