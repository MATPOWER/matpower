Change history for MOST
=======================


Version 1.0.1 - *Oct 30, 2016*
------------------------------

#### 10/30/18
  - Release 1.0.1.

#### 10/26/18
  - **INCOMPATIBLE CHANGE**: Failure of the optimization no longer
    halts execution and jumps to the debugger.
  - Add `success` flag to `md.results` output MOST Data struct to
    indicate success or failure of optimization.

#### 3/19/18
  - Fix bugs in `plot_uc_data()` resulting in incorrect legends.

#### 3/7/18
  - Replace `clock()`/`etime()` with `tic()`/`toc()` for timing.

#### 9/27/17
  - Fix bug #1 in `loadmd()` where profiles that modify xGenData or
    StorageData resulted in a fatal error.

#### 9/12/17
  - Update for `@opt_model` API cleanup.

#### 9/1/17
  - Use MATPOWER's new quadratic costs in `@opt_model` in place of the
    legacy cost model.

#### 8/7/17
  - Switch to OOP notation everywhere for `opt_model` object,
    e.g. `om.method()`.
  - Use `om.init_indexed_name()` instead of deprecated form of calls to
    `add_vars()`, `add_constraints()` or `add_costs()`.
  - Use `om.add_lin_constraint()` in place of deprecated `add_constraints()`.

#### 5/25/17
  - Fix dimension of `RampWear` cost indexing if `mdi.OpenEnded` is true.
  - Add missing constant term to objective function value reported
    by `most_summary`.

#### 1/26/17
  - Add MOST User's Manual to `docs` and sources to `docs/src`.

#### 12/21/16
  - Add Travis-CI integration. *Thanks to Richard Lincoln.*


Version 1.0 - *Dec 16, 2016*
----------------------------

#### 12/16/16
  - Released 1.0.
  - Moved development to GitHub: <https://github.com/MATPOWER/most>.

#### 3/1/16
  - Put checks in `loadmd()` and `most()` to require internal bus ordering
    for `mpc`, plus explicit notes in documentation.

#### 2/26/16
  - Rename example files in `most/t` to start with `ex_` instead of `eg_`.

#### 2/25/16
  - Add MATPOWER option `most.fixed_res` with default of -1 (depends
    on presence of `md.FixedReserves`) to control `md.IncludeFixedReserves`.
  - Move fixed zonal reserve output fields (`R`, `prc`, `mu.l`, `m.u`,
    `m.Pmax`, `totalcost`) from `mdo.FixedReserves(t,j,k)` to
    `mdo.flow(t,j,k).mpc.reserves`.

#### 2/24/16
  - Add `most_summary()` function to summarize and print some summary
    results. Moved from some of the test files into its own public
    function. (Still very incomplete).
  - For consistency, rename MOST data struct variables everywhere to
    `mdi` (input) from `Istr`, `md_in`, and `mdin`, and to `mdo` (output)
    from `Ostr`, `md_out`, and `mdout`.
  - Rename `md.idx` fields for dimensions of dynamic system constraints.
    - `nyo`  -->  `nyds`
    - `nzd`  -->  `nzds`
    - `nyt`  -->  `ntds`

#### 2/22/16
  - (Tentatively) modify `md_init()` to initialize only what is needed to
    run `loadmd()` and `most()`.

#### 1/26/16
  - Renamed `mops` to `most`, and re-wrote history below accordingly.

#### 12/18/15
  - Fixed fatal crash triggered by failed solve with `verbose` option off.

#### 12/9/15
  - Added `ExpectedTerminalStorageMax`,` ExpectedTerminalStorageMin` to
    `md.Storage` and `StorageData` structs. If present
    `ExpectedTerminalStorageAim` now simply overwrites both.

#### 10/27/15
  - Renamed `apply_contingency()` to `apply_changes()` and moved from sopf
    to matpower.

#### 7/17/15
  - Fixed bug preventing proper printing of exit flag on failed solve.

#### 7/10/15
  - Added `mpopt` as optional second argument to `most()`.

#### 7/9/15
  - Removed all of the old indexing fields from `md_init()`.
  - Added `most` to `have_fcn()`.
  - Added `most` options to `mpoption()`, added `mpoption_info_most()`.

#### 7/2/15
  - Renamed `loadmpsd()` to `loadmd()`, `mpsd_init()` to `md_init()`,
    `mpsd` to `md`.

#### 7/1/15
  - Added `plot_uc()` function for plotting unit commitment schedules.

#### 6/19/15
  - Moved lots of `mpsopf` related files to `matpower/dist/most` temporarily.
  - Commented out call to `oldidx()`.

#### 6/12/15
  - In `loadmpsd()` `xgd.CommitKey` must be non-empty (in addition to just
    being present) in order to process unit commitment fields.

#### 6/10/15
  - Added `fixed_gencost` field to flow-specific `mpc` fields to save
    the fixed cost portion that is removed from `gencost`, to allow
    for computation of full flow-specific generator costs from solution.

#### 6/9/15
  - Modified `filter_ramp_transitions()` to multiply conditional probability
    of transition by conditional probability of being in source state
    (assuming we've made it to that period) before applying cutoff
    threshold. This will typically cut off more transitions for the same
    threshold value.
  - Minor fixes, updates to `plot_gen()`.

#### 5/21/15
  - Forked development from `mpsopfl_fixed_res()`, which is currently
    identical except for function names in error messages.
