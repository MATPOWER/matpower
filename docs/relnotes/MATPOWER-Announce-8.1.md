What's New in MATPOWER 8.1
--------------------------

#### Released Jul 12, 2025

Below are some of the highlights of the changes since version 8.0 of
MATPOWER. See the [full release notes][1] and the [`CHANGES.md`][2]
file for more details. For release notes for previous versions, see
Appendix H of the [MATPOWER User's Manual][3].


#### New Features:

- Three-Phase Proof-of-Concept Enhancements
  - New utility converts a standard single-phase MATPOWER case to an equivalent balanced three-phase case.
  - New prototype three-phase shunt model and off-nominal taps in transformer model.
  - Support for three-phase prototype data in `savecase()`.
- [MP-Opt-Model][4] 5.0 includes new, redesigned classes for building and solving mathematical programming/optimization models, support for quadratic constraints and quadratically constrained quadratic programming (QCQP) models and solvers, support for the open-source [HiGHS][5] solver for LP, QP and MILP problems, and new `relax_integer` option for solving the LP/QP relaxation of MILP/MIQP model.
- MOST Pro 1.4.1, available by contacting [info@matpower.org][6], adds support for MATPOWER DC lines as described in Section 7.6.3 of the [MATPOWER User's Manual][3]. _Please note:_ MATPOWER 8.1 includes MOST 1.3.1, and MOST Pro 1.4.1 is a paid upgrade. 
- New functions:
  - `mp.case_utils.convert_1p_to_3p()` -- Converts a single-phase case to a three-phase case.
- New options:
  - `highs.opts` overrides default options for the [HiGHS][5] solver.


#### New Case Files:

- Two new cases.
  - `case1197` -- 1197-bus radial distribution system case.
  - `case59` -- 59 bus, 14 generator Australian network case.


#### New Documentation:
Three live scripts illustrate the use of new features.
- single-phase to three-phase conversion capabilities
- using new `mp.opt_model` class to build and solve an optimization model
- new quadratic constraint features and QCQP features


#### Other Improvements:

- Update versions of included packages:
  - MP-Test 8.1
  - MIPS 1.5.2
  - MP-Opt-Model 5.0
  - MOST 1.3.1
- Numerous bug fixes.

*Thanks to Wilson González Vanegas for numerous contributions to this release.*

[1]: https://github.com/MATPOWER/matpower/blob/master/docs/relnotes/MATPOWER-Release-Notes-8.1.md
[2]: https://github.com/MATPOWER/matpower/blob/master/CHANGES.md
[3]: https://github.com/MATPOWER/matpower/blob/master/docs/MATPOWER-manual.pdf

[4]: https://github.com/MATPOWER/mp-opt-model
[5]: https://highs.dev
[6]: mailto:info@matpower.org?subject=MOST%20Pro&body=Please%20send%20me%20information%20on%20obtaining%20MOST%20Pro
