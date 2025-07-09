What's New in MOST 1.3.1
------------------------

#### Released Jul 7, 2025

Below is a summary of the changes since version 1.3 of MOST. See the
[`CHANGES.md`][1] file for all the gory details. For release notes for
previous versions, see Appendix B of the [MOST User's Manual][2].

#### Changes:
  - `most_summary()` now skips display of non-existent contingencies.
  - Include [HiGHS](https://highs.dev) solver (via
    [HiGHSMEX](https://github.com/savyasachi/HiGHSMEX) interface), if
    available, in unit commitment tests.

#### Bugs Fixed:
  - Fix [issue #45][3] where `most()` does not properly handle cases with
    contingencies defined only in some periods/scenarios.
    *Thanks to Stefano Nicolin.*
  - Fix issue with `most_summary()` when ramp results are missing.
  - Tweak tests to work around bug in HiGHS-based `linprog` and
    `intlinprog` in Optimization Toolbox R2024a and R2024b.
  - Fix issue caused by tiny non-zero values for commitment variables.
    Don't count on MP-Opt-Model's `miqps_<solver>()` functions to round
    integer variable solutions.

#### PRO Version
  - There is now a PRO version of MOST that adds support for MATPOWER DC
    lines as described in Section 7.6.3 of the [MATPOWER User's Manual][4].
    Please contact [info@matpower.org](mailto:info@matpower.org?subject=MOST%20Pro&body=Please%20send%20me%20information%20on%20obtaining%20MOST%20Pro.)
    for information on obtaining this version.


[1]: ../../CHANGES.md
[2]: ../MOST-manual.pdf
[3]: https://github.com/MATPOWER/most/issues/45
[4]: https://matpower.org/docs/MATPOWER-manual-8.1.pdf
