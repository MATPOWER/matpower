What's New in MOST 1.3
----------------------

#### Released May 10, 2024

Below is a summary of the changes since version 1.2 of MOST. See the
[`CHANGES.md`][1] file for all the gory details. For release notes for
previous versions, see Appendix B of the [MOST User's Manual][2].

#### New Features:
  - New Sphinx-based [Reference documentation][3].

#### Changes:
  - Reduce memory requirements for long horizon cases with storage by
    forming/storing transposes of matrices for storage constraints.
    *Requires [MP-Opt-Model][4] 4.2 or later.*
  - Speed up building unit commitment (min up/down time) constraints.
    Improvement can be quite substantial on large problems.

#### Bugs Fixed:
  - Fix [issue #37][5] which caused a fatal error in storage input
    checks with multiple storage units under some circumstances.
    *Thanks to Keir Steegstra.*
  - Fix [issue #39][6] in which the value of `mdi.Delta_T`, the number
    of hours represented by each period, was not being accounted for in
    most of the terms in the objective function.
    *Thanks to Stefano Nicolin.*

#### Incompatible Changes:
  - Remove extra column in `mdo.results.ExpectedRampCost` and ignore for
    single period.


[1]: ../../CHANGES.md
[2]: ../MOST-manual.pdf
[3]: https://matpower.org/doc/most/
[4]: https://github.com/MATPOWER/mp-opt-model/
[5]: https://github.com/MATPOWER/most/issues/37
[6]: https://github.com/MATPOWER/most/issues/39
