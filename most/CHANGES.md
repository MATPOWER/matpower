Change history for MOST
=======================


Version 1.3 - *May 10, 2024*
----------------------------

#### 5/10/24
  - Release 1.3.

#### 4/4/24
  - Tweak tests to avoid warnings and presolve bug with new HiGHS-based
    `linprog` and `intlinprog` in Optimization Toolbox R2024a.

#### 3/21/24
  - Add Sphinx-based Reference documentation.

#### 11/8/23
  - Speed up building unit commitment (min up/down time) constraints.
    Improvement can be quite substantial on large problems.

#### 11/6/23
  - Reduce memory requirements for long horizon cases with storage by
    forming/storing transposes of matrices for storage constraints.
    _Requires [MP-Opt-Model][11] version > 4.1._

#### 10/25/23
  - Fix [issue #39][10] in which the value of `mdi.Delta_T`, the number of
    hours represented by each period, was not being accounted for in most
    of the terms in the objective function.
    *Thanks to Stefano Nicolin.*

#### 10/4/23
  - Fix [issue #37][9] which caused a fatal error in storage input checks
    with multiple storage units under some circumstances.
    *Thanks to Keir Steegstra.*

#### 2/3/23
  - Remove extra column in mdo.results.ExpectedRampCost and ignore for
    single period.


Version 1.2 - *Dec 13, 2022*
----------------------------

#### 12/13/22
  - Release 1.2.

#### 11/10/22
- Update `addgen2mpc()` to expand fixed reserve parameters to keep sizes
  compatible.

#### 9/29/22
- Silence near singular matrix warnings in some tests that began with
  MATLAB R2022b.

#### 8/30/22
  - Fix [issue #29][8], where a typo caused a check on
    `md.UC.MinDown` >= 1 to be skipped.
    *Thanks to Talha Iqbal.*

#### 6/15/22
  - Add TLMP (temporal locational marginal price) for storage units
    based on work by Chen, Tong in [[Chen2022]][7] and returned in
    `mdo.results.StorageTLMPc`, `mdo.results.StorageTLMPd`, 
    `mdo.results.CondStorageTLMPc`, `mdo.results.CondStorageTLMPd`.
    See Table 5-13 in the User's Manual.
  - Add tests for TLMP in `t_most_tlmp` based on toy examples from Cong Chen,
    including both ramping and storage.

#### 6/14/22
  - For deterministic cases with storage where `ForceCyclicStorage` is 0,
    ensure that initial storage bounds are equal to initial storage and
    output a warning if they are modified. Fix deterministic UC tests
    where this was causing results to change depending on value of `rho`.

#### 5/12/22
  - Add calculation of expected TLMP (temporal locational marginal price)
    based on work by Guo, Chen, Tong in [[Guo2021]][5] and [[Chen2021]][6]
    and returned in `mdo.results.GenTLMP` and `mdo.results.CondGenTLMP`.
    See Table 5-13 in the User's Manual.

#### 5/4/22
  - Ramping reserves and constraints are now included for the transition
    from the initial state into period 1, except for single-period problems.
  - **INCOMPATIBLE CHANGE**: Modified definition of ramping reserves for
    period _t_ (and all corresponding input and output parameters) to refer
    to the transition from _t-1_ to _t_, not _t_ to _t+1_. This means that
    the ramping reserves for the transition into the first period are now
    optimization variables and the corresponding constraints are explicit.
    This is for multiperiod problems only. Ramping reserves and contraints
    are explicitly excluded for single-period problems.  
    _Note:_ This change also corrects an error in (4.11) in the manual. The
    superscript _t_ on gamma is now correct. Previously it should have been
    _t+1_, as it was in the code.

#### 4/22/22
  - Fix tests that were failing under Octave 7.x.

#### 5/24/21
  - Fix bug in `plot_uc()` that prevented it from working properly in a
    subplot, such as in `t_most_uc()`.
    *Thanks to Lim Han.*


Version 1.1 - *Oct 8, 2020*
---------------------------

*Requires MATPOWER with 7.1 or later (for MP-Opt-Model 3.0 or later).*

#### 10/8/20
  - Release 1.1.

#### 9/9/20
  - Use `opt_model.get_soln()` to extract variable and shadow price
    results, rather than doing the indexing manually.

#### 3/19/20
  - Convert to using `opt_model.solve()` method rather than calling
    `miqps_matpower()` or `qps_matpower()` directly.
  - **INCOMPATIBLE CHANGE**: Update objective function value returned in
    `mdo.QP.f` to include the previously missing constant term.

#### 2/27/20
  - Fix [issue #16][4], where the `om` field of the output MOST data
    struct (`mdo`) was a handle to the same object as as the `om`
    field of the input MOST data struct (`mdi`), meaning that changing
    one would modify the other.
    *Thanks to Baraa Mohandes.*

#### 8/27/19
  - Update `most_summary` to include sections for fixed loads and
    storage expected stored energy.
    *Thanks to Baraa Mohandes.*
  - Move assembly of constraints and variable bounds inside the
    `build_model` conditional.
    *Thanks to Baraa Mohandes.*

#### 8/21/19
  - Fix [bug #6][3] where building a model without solving it, or
    solving a previously built model resulted in a fatal error.
    *Thanks to Baraa Mohandes.*

#### 8/20/19
  - Fix [bug #11][2] where storage constraints were not correct for
    t=1 and `rho ~= 1`. *Thanks to Baraa Mohandes.*


Version 1.0.2 - *Jun 20, 2019*
------------------------------

#### 6/20/19
  - Release 1.0.2.
  - Add `CITATION` file.
  - Other miscellaneous documentation updates, e.g. MATPOWER website
    links updated to https://matpower.org, separate references for
    MATPOWER software and User's Manual, with DOIs.

#### 11/20/18
  - Fix selection of default solver for `t_most_w_ds` and add option
    to specify solver directly in optional input arg.


Version 1.0.1 - *Oct 30, 2018*
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
  - Fix [bug #1][1] in `loadmd()` where profiles that modify xGenData
    or StorageData resulted in a fatal error.

#### 9/12/17
  - Update for `opt_model` API cleanup.

#### 9/1/17
  - Use MATPOWER's new quadratic costs in `opt_model` in place of the
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
  - no changes from v1.0b2

#### 11/17/16
  - Silence some warnings during tests on Octave 4.2.


Version 1.0b2 - *Nov 1, 2016*
-----------------------------

#### 11/1/16
  - Released 1.0b2.

#### 10/27/16
  - Fixed some MOSEK related issues in tests and tutorial examples.


Version 1.0 - *Jun 1, 2016*
---------------------------

#### 6/1/16
  - Released 1.0b1.

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

---

[1]: https://github.com/MATPOWER/most/issues/1
[2]: https://github.com/MATPOWER/most/issues/11
[3]: https://github.com/MATPOWER/most/issues/6
[4]: https://github.com/MATPOWER/most/issues/16
[5]: https://doi.org/10.1109/TPWRS.2021.3055730
[6]: https://doi.org/10.1109/TPWRS.2020.3045162
[7]: https://arxiv.org/abs/2204.08140
[8]: https://github.com/MATPOWER/most/issues/29
[9]: https://github.com/MATPOWER/most/issues/37
[10]: https://github.com/MATPOWER/most/issues/39
[11]: https://github.com/MATPOWER/mp-opt-model
