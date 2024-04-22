What's New in MOST 1.1
----------------------

#### Released Oct 8, 2020

Below is a summary of the changes since version 1.0.2 of MOST. See the
[`CHANGES.md`][1] file for all the gory details. For release notes for
previous versions, see Appendix B of the [MOST User's Manual][2].

#### Changes:
  - Requires [MATPOWER][4] 7.1 or later.
  - Output of `most_summary()` includes sections for fixed loads and for
    expected stored energy for storage units.
  - Relies on [MP-Opt-Model][3] 3.0, which can be found at
    https://github.com/MATPOWER/mp-opt-model and is included in
    [MATPOWER][4] 7.1.
    - Significant performance improvement for some problems when constructing
    sparse matrices for linear constraints or quadratic costs (e.g. during
    problem setup).
    *Thanks to Daniel Muldrew.*
    - Uses the `opt_model.solve()` method rather than calling
      `miqps_matpower()` or `qps_matpower()` directly.
    - Uses the `opt_model.get_soln()` method to extract variable and shadow
      price results, rather than doing the indexing manually.

#### Bugs Fixed:
  - Fix bug [#6][5] where building a model without solving it, or solving a
    previously built model resulted in a fatal error.
    *Thanks to Baraa Mohandes.*
  - Fix bug [#11][6] where storage constraints were not correct for
    `t=1` and `rho` not equal to 1. *Thanks to Baraa Mohandes.*
  - Fix issue [#16][7], where the `om` field of the output MOST data struct
    (`mdo`) was a handle to the same object as as the `om` field of the
    input MOST data struct (`mdi`), meaning that changing one would modify
    the other. *Thanks to Baraa Mohandes.*

#### Incompatible Changes:
  - Objective function value returned in `mdo.QP.f` updated to include the
    previously missing constant term.


[1]: ../../CHANGES.md
[2]: ../MOST-manual.pdf
[3]: https://github.com/MATPOWER/mp-opt-model
[4]: https://github.com/MATPOWER/matpower
[5]: https://github.com/MATPOWER/most/issues/6
[6]: https://github.com/MATPOWER/most/issues/11
[7]: https://github.com/MATPOWER/most/issues/16
