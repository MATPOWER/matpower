What's New in MATPOWER 7.1
--------------------------

#### Released Oct 8, 2020

Below are some of the highlights of the changes since version 7.0 of
MATPOWER. See the [full release notes][1] and the [`CHANGES.md`][2]
file for more details. For release notes for previous versions, see
Appendix H of the [MATPOWER User's Manual][3].

#### New Features:
- Core optimization model and solver interfaces released as separate
  package (included), [MP-Opt-Model][4] 3.0, bringing many enhancements:
  - Support for new QP solver [OSQP][5] for DC OPF
  - New unified interfaces for nonlinear programming and nonlinear equation
    solvers
  - Many new features of the `opt_model` class, including a `solve()`
    method that selects and calls the appropriate solver based on the
    automatically detected problem type, and methods for extracting
    specific variables, costs, constraint values and shadow prices, etc.
    from a solved model
  - Performance improvements
- [MP-Test][6] 7.1, with new modular, extensible `have_feature()` function
  for detecting optional functionality
- New options and speed improvements for computing shift factors.
- Numerous new functions and program options.

#### New Case Files:
- Nineteen new distribution system cases \[[ref 1][7], [ref 2][8]\].
  *Thanks to Houssem Bouchekara, et. al.*

#### New Documentation:
- New [MP-Opt-Model User's Manual][9], included in `mp-opt-model/docs`.

#### Other Improvements:
- New implementation of Newton power flow for cartesian voltage
  representations.
- Improved robustness of Newton power flow with hybrid voltage updates.
- Significant performance improvement for CPLEX on small DC OPF problems.
- Update versions of included packages:
  - MIPS 1.4.
  - MOST 1.1.
  - MP-Test 7.1.
- Numerous bug fixes.

#### Incompatible Changes:
- Update `case18`, `case22`, `case69`, `case85` and `case141` to more closely
  match data from original papers, thanks in part to case files contributed
  by Houssem Bouchekara, et al. Solutions for updated cases may not match
  exactly. See help text in case files for details.


[1]: MATPOWER-Release-Notes-7.1.md
[2]: ../../CHANGES.md
[3]: ../MATPOWER-manual.pdf
[4]: https://github.com/MATPOWER/mp-opt-model
[5]: https://osqp.org
[6]: https://github.com/MATPOWER/mptest
[7]: https://doi.org/10.18799/24056537/2019/3/196
[8]: https://doi.org/10.36227/techrxiv.12578648.v1
[9]: https://matpower.org/docs/MP-Opt-Model-manual-3.0.pdf
