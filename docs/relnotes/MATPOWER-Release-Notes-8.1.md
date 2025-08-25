What's New in MATPOWER 8.1
--------------------------

#### Released Jul 12, 2025

Below is a summary of the changes since version 8.0 of MATPOWER. See the
[`CHANGES.md`][1] file for all the gory details. For release notes for
previous versions, see Appendix H of the [MATPOWER User's Manual][2].


#### New Features:
- Three-Phase Proof-of-Concept Enhancements
  - New utility, `mp.case_utils.convert_1p_to_3p()`, converts a standard single-phase MATPOWER case to an equivalent balanced three-phase case.
    *Thanks to Wilson González Vanegas.*
  - New off-nominal tap ratio parameter in prototype three-phase transformer model.
    *Thanks to Wilson González Vanegas.*
  - New prototype three-phase shunt model.
    *Thanks to Wilson González Vanegas.*
  - Support for three-phase prototype data in `savecase()`.
- [MP-Opt-Model][3] 5.0 includes:
  - New, redesigned classes for building and solving mathematical programming/optimization models. The new classes `mp.opt_model` and `mp.set_manager` and subclasses replace refactored `opt_model` and `mp_idx_manager` which are deprecated, but retained for backward compatibility.
  - Support for quadratic constraints and quadratically constrained quadratic programming (QCQP) models and solvers.
    *Thanks to Wilson González Vanegas.*
  - Support for the open-source [HiGHS][4] solver for LP, QP and MILP problems based on the [HiGHSMEX][5] interface.
  - New `relax_integer` option to facilitate easily solving the LP or QP relaxation of a mixed integer LP or QP model.
  - Support for [Artelys Knitro][6] solver for LP and QP problems.
    *Thanks to Wilson González Vanegas.*
  - Support for [Artelys Knitro][6] 15.x.
  - Support for conversion between objects and structs to facilitate workarounds for Octave's inability to save and load classdef objects.
  For details, see the [MP-Opt-Model 5.0 release notes][7].
- MOST Pro 1.4.1, available by contacting [info@matpower.org][8], adds support for MATPOWER DC lines as described in Section 7.6.3 of the [MATPOWER User's Manual][9]. _Please note:_ MATPOWER 8.1 includes MOST 1.3.1, and MOST Pro 1.4.1 is a paid upgrade. 
- [MP-Test][10] 8.1 adds functions to assist with code for debugging and enhances handling of `Inf` values in test results. For details, see the [MP-Test 8.1 release notes][11].
- New functions:
  - `mp.case_utils.convert_1p_to_3p` -- Converts a single-phase case to a three-phase case.
  - `mp.load_dm` -- Loads a data model object.
  - `save2psse_rop` -- Saves data from a MATPOWER case to a PSS/E ROP (Raw Operating Point) file. *Thanks to Irabiel Romero.*
- New options:
  - `highs.opts` overrides default options for the [HiGHS][4] solver.


#### New Case Files:
- One new distribution system case. *Thanks to Paul S. Moses.*
  - `case1197` -- 1197-bus radial distribution system case.
- One new transmmission system case. *Thanks to Mory Najafi.*
  - `case59` -- 59 bus, 14 generator Australian network case.


#### New Documentation:
Three live scripts illustrate the use of new features.
- `convert_1p_to_3p_ex1.mlx` (in `examples`) illustrates the use of the new single-phase to three-phase conversion capabilities. [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=MATPOWER/matpower&project=matpower.prj&file=examples/convert_1p_to_3p_ex1.mlx)
- `milp_example1.mlx` (in `mp-opt-model/examples`) illustrates the use of MP-Opt-Model and the new `mp.opt_model` class to build and solve an optimization (MILP) model. [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=MATPOWER/matpower&project=matpower.prj&file=mp-opt-model/examples/milp_example1.mlx)
- `qcqp_example1.mlx` (in `mp-opt-model/examples`) illustrates the new quadratic constraint features and two methods of building and solving a quadratically-constrained quadratic programming (QCQP) model. [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=MATPOWER/matpower&project=matpower.prj&file=mp-opt-model/examples/qcqp_example1.mlx)


#### Other Improvements:
- Add shunt loss columns (`psh_fr`, `qsh_fr`, `psh_to`, `qsh_to`) to branch data model. In pretty-printed output, separate branch loss summary into series and shunt losses and add series loss columns to branch detail output.
*Thanks to Wilson González Vanegas.*
- The `save2psse` function now creates distinct generator IDs for cases with multiple generators at a single bus. *Thanks to Irabiel Romero.*
- [MIPS][12] 1.5.2 adds a feature detection function for \mips{}. For details, see the [MIPS 1.5.2 release notes][13].
- [MOST][14] 1.3.1 fixes bugs in `most_summary()` and some tests, and includes unit commitment tests for [HiGHS][4] solver. For details, see the [MOST 1.3.1 release notes][15].
- Include [HiGHS][4] solver, if available, in DC OPF tests.
- Add `rebuild()` method to `mp.data_model`, used to update the data model internals after making changes to data in the data tables of the elements.
- Make optional the `mpopt` argument to `mp.task.run()`.


#### Bugs Fixed:
- Fix typos that cause fatal errors with [GNU Octave][16] 10.x.
- Fix typo in units of min and max reactive power shadow price in summary section of legacy OPF output.
-  Optional `'soln_fname'` input to `run_mp()`, hence also to `run_pf()`, `run_cpf()`, and `run_opf()`, was being ignored.
- Fix bug #256 to make `opf.use_vg` option compatible with generators at PQ buses.
- Build load and shunt connectivity using `bus` column in main table (i.e. bus numbers), not `source_uid` column (i.e. row indices in original `bus` table).


#### Incompatible Changes:
- See the [MP-Opt-Model 5.0 release notes][7] for incompatible changes related to [MP-Opt-Model][3].


[1]: https://github.com/MATPOWER/matpower/blob/master/CHANGES.md
[2]: https://github.com/MATPOWER/matpower/blob/master/docs/MATPOWER-manual.pdf
[3]: https://github.com/MATPOWER/mp-opt-model
[4]: https://highs.dev
[5]: https://github.com/savyasachi/HiGHSMEX
[6]: https://www.artelys.com/solvers/knitro/
[7]: https://github.com/MATPOWER/mp-opt-model/blob/master/docs/relnotes/MP-Opt-Model-Release-Notes-5.0.md
[8]: mailto:info@matpower.org?subject=MOST%20Pro&body=Please%20send%20me%20information%20on%20obtaining%20MOST%20Pro
[9]: https://matpower.org/docs/MATPOWER-manual-8.1-dev.pdf
[10]: https://github.com/MATPOWER/mptest
[11]: https://github.com/MATPOWER/mptest/blob/master/docs/relnotes/MP-Test-Release-Notes-8.1.md
[12]: https://github.com/MATPOWER/mips
[13]: https://github.com/MATPOWER/mips/blob/master/docs/relnotes/MIPS-Release-Notes-1.5.2.md
[14]: https://github.com/MATPOWER/most
[15]: https://github.com/MATPOWER/most/blob/master/docs/relnotes/MOST-Release-Notes-1.3.1.md
[16]: https://www.octave.org