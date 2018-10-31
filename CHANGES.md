Change history for MATPOWER
===========================

MATPOWER release notes, found in [`docs/relnotes`](docs/relnotes),
contain summaries of the primary changes in each versioned release.

For change history for [MP-Test][1], see [mptest/CHANGES.md](mptest/CHANGES.md).  
For change history for [MIPS][2], see [mips/CHANGES.md](mips/CHANGES.md).  
For change history for [MOST][3], see [most/CHANGES.md](most/CHANGES.md).


Version 7.0b1 - *Oct 31, 2018*
------------------------------

#### 10/31/18
  - Released 7.0b1.
  - Update versions of included packages
    - MIPS 1.3
    - MOST 1.0.1
    - MP-Test 7.0b1

#### 10/25/18
  - Add another purely synthetic case from the ACTIVSg team (ASU, Cornell,
    Texas A&M, U of IL, and VCU - Synthetic grids), resulting from work
    supported by the ARPA-E GRID DATA program.
    *Thanks to Adam Birchfield and the ACTIVSg team.*
    - `case_SyntheticUSA` (82,000-bus continental USA synthetic model,
      aggregation of `case_ACTIVSg2000`, `case_ACTIVSg10k`,
      `case_ACTIVSg70k`, connected by 9 DC lines)

#### 10/16/18
  - **INCOMPATIBLE CHANGE:** Correct signs of phase shifter angles in
    Polish system cases, since they were based on the old sign convention
    used by MATPOWER prior to v3.2 (see change on 6/21/07). Affects the
    following cases:
    - `case2383wp`
    - `case2736sp`
    - `case2737sop`
    - `case2746wop`
    - `case2746wp`
    - `case3375wp`
    *Thanks to Mikhail Khokhlov and Dr. Artjoms Obusevs for reporting.*
  - Fix `runpf()` handling of case where individual power flow fails
    during Q limit enforcement.

#### 10/9/18
  - Add option `opf.softlims.default` for use with enhanced
    `toggle_softlims()` to determine whether or not to include soft limits
    on constraints whose parameters are not specified explicitly in the
    `mpc.softlims` struct.
  - Update `toggle_softlims()` to implement soft limits for all OPF inequality
    constraints, i.e. bus voltage bounds, generator active & reactive bounds,
    branch flow and branch angle difference limits.
    *Thanks to Eran Schweitzer.*
  - Add optional `mpopt` argument to calls to `ext2int()` and `int2ext()`.
  - **INCOMPATIBLE CHANGE:** Add `mpopt` to input args for OPF `ext2int`
    and `int2ext` callbacks.
  - **INCOMPATIBLE CHANGE:** Turning soft limits on without specifying any
    parameters explicitly in `mpc.softlims` now implements soft limits for
    all constraints, by default, not just branch flow limits. And the
    format of the input parameters in `mpc.softlims` has changed. See
    `help toggle_softlims` for the details.

#### 9/10/18
  - Add support for PARDISO v6.x.

#### 8/30/18
  - Update to Aug 9, 2018 version of `case_ACTIVSg2000`.

#### 6/19/18
  - Correct (swap) names of fields `e2i` and `i2e` used internally
    (in `mpc.order.gen`) to map between external and internal indexing
    of generators.
  - Improve error messages and documentation related to specifying
    dispatchable loads and their power factor constraints.

#### 5/31/18
  - Add another purely synthetic case from the ACTIVSg team (ASU, Cornell,
    Texas A&M, U of IL, and VCU - Synthetic grids), resulting from work
    supported by the ARPA-E GRID DATA program.
    *Thanks to Adam Birchfield and the ACTIVSg team.*
    - `case_ACTIVSg70k` (70,000-bus Eastern US synthetic model)

#### 5/19/18
  - Add `loadshed()` function to compute MW curtailments of dispatchable
    loads.

#### 5/18/18
  - Add `lng` for "Liquefied Natural Gas" to `genfuels()`.
  - Under Octave, disable the use of explicit LU decomposition with AMD
    reordering and the 3 output argument form of LU for Newton power
    flow on larger systems. In Octave 4.4 it generates warnings about
    potential LU failures **and** slows it down.

#### 5/10/18
  - Fix bug in order of lambdas used in `opf_branch_ang_hess()` and
    `opf_vlim_hess()` that resulted in wrong sign in Hessian terms.
  - Fix bug #36 where Hessian structure for Ipopt and Knitro were
    incorrect. Re-enabled skipped tests that were previously failing.

#### 4/19/18
  - Add AC OPF tests for cases with ref bus ~= 1 and ref angle ~= 0.
  - Fix fatal error in `opf_vref_hess()` for cases where ref bus ~= 1.

#### 4/12/18
  - Add to `nested_struct_copy()` ability to copy fields that are struct
    arrays.

#### 4/5/18
  - **INCOMPATIBLE CHANGE:** Swap the order of the output arguments of
    `dSbus_dV()` for polar coordinate voltages (angle before magnitude)
    for consistency.

#### 4/3/18
  - Add three new variants of the standard AC OPF formulation, for a total
    of four, including both nodal power and current balance constraints and
    both polar and cartesian representations of voltage. See the new
    `opf.current_balance` and `opf.v_cartesian` options.
    *Thanks to Baljinnyam Sereeter.*
  - Update existing functions that compute derivatives with respect to
    voltage to handle cartesian coordinate voltages in addition to polar.
    *Thanks to Baljinnyam Sereeter.*
  - Update existing functions that compute derivatives with respect to
    voltage to handle cartesian coordinate voltages in addition to polar.
    *Thanks to Baljinnyam Sereeter.*
  - Add functions to compute derivatives of current balance constraint:
    *Thanks to Baljinnyam Sereeter.*
    - `dImis_dV`
    - `d2Imis_dV2`
    - `d2Imis_dVdSg`
  - Add functions to implement constraints for new AC OPF variants.
    *Thanks to Baljinnyam Sereeter.*
    - `opf_branch_ang_fcn`, `opf_branch_ang_hess`
    - `opf_current_balance_fcn`, `opf_current_balance_hess`
    - `opf_veq_fcn`, `opf_veq_hess`
    - `opf_vlim_fcn`, `opf_vlim_hess`
    - `opf_vref_fcn`, `opf_vref_hess`
  - Add two new Technical Notes (TN3, TN4) and updated revision of TN2.
    *Thanks to Baljinnyam Sereeter.*
    - [MATPOWER Technical Note 3](http://www.pserc.cornell.edu/matpower/TN3-More-OPF-Derivatives.pdf)
    - [MATPOWER Technical Note 4](http://www.pserc.cornell.edu/matpower/TN4-OPF-Derivatives-Cartesian.pdf)
  - Add `d2Abr_dV2()` to replace the (now deprecated) separate functions
    `d2AIbr_dV2()` and `d2ASbr_dV2()`.

#### 3/7/18
  - Replace clock()/etime() with tic()/toc() for timing.

#### 12/22/17
  - Add another purely synthetic case from the ACTIVSg team (ASU, Cornell,
    Texas A&M, U of IL, and VCU - Synthetic grids), resulting from work
    supported by the ARPA-E GRID DATA program.
    *Thanks to Adam Birchfield and the ACTIVSg team.*
    - `case_ACTIVSg25k` (25,000-bus US Northeast/Mid-Atlantic synthetic
       model)

#### 12/6/17
  - Fix bug #21 where a continuation power flow that failed the first
    corrector step would produce a fatal error.
    *Thanks to Elis Nycander.*
  - Fix bug #23 where the continuation power flow could switch
    directions unexpectedly when the operating point switched from
    stable to unstable manifold or vice-versa after hitting a limit.
    *Thanks to Elis Nycander and Shrirang Abhyankar.*
  - Fix bug #26 where, in a continuation power flow, a reactive limit
    at a bus could be detected in error if multiple generators at the bus
    had reactive ranges of very different sizes.
    *Thanks to Elis Nycander and Shrirang Abhyankar.*
  - Remove numerical proxies for `Inf` from case files. Some cases used
    9999 and/or 99999 as numerical proxies to indicate the absence of
    limits on generator `QMIN`, `QMAX` (replaced with `-Inf`, `Inf`) or
    branch `RATE_A`, `RATE_B`, `RATE_C` (replaced with 0).

#### 11/29/17
  - Add new option `knitro.maxit` to set maximum number of iterations
    for AC OPF solver using Knitro.

#### 11/22/17
  - Add new option `opf.start` to replace deprecated `opf.init_from_mpc`
    and add a new possibility to automatically run a power flow to
    initialize the starting state for the OPF.

#### 11/17/17
  - Improve handling of allocation of reactive power among multiple
    generators at a bus in power flow solution. Use equal violation at
    buses where total Qmin == total Qmax. Eliminate NaNs in case of
    infinite Q ranges.

#### 11/14/17
  - Add MATPOWER logo to User's Manual.

#### 11/7/17
  - Remove deprecated non-fatal error mechanism in `loadcase()`.
  - Update `feval_w_path()` to allow empty string for path, equivalent to
    calling `feval()` directly.
  - Modify `loadcase()` to use `feval_w_path()`.
  - Make `install_matpower()` check for minimum required MATLAB/Octave version.
  - __IMPORTANT NOTE__: For MATLAB users, the minimum requirement for
    MATPOWER is now MATLAB 7.3.0 (R2006b). _(No access to older versions
    for testing)_.

#### 11/3/17
  - Add `case_RTS_GMLC` from [here](https://github.com/GridMod/RTS-GMLC).

#### 10/31/17
  - Update `case9` with complete reference to source data and updated
    generator voltages, initial reactive injections and slack active
    injection to more closely match original source.
  - Make AC OPF solver always default to MIPS, even if TSPOPF is installed.
    i.e. `opf.ac.solver = 'DEFAULT'` is now identical to
    `opf.ac.solver = 'MIPS'`.
  - Add [E4ST](http://e4st.com/) to `have_fcn()` and `mpver` output.

#### 10/27/17
  - Add contingencies and scenarios for ACTIVSg cases, contingencies for
    all four, and one year of hourly zonal or area load scenarios for 200
    and 2000 bus cases.

#### 10/11/17
  - Add `t_opf_default()` to test AC OPF with `opf.ac.solver` set to `DEFAULT`.
  - Fix bug in setup of OPF (introduced since v6.0), triggered by running
    an AC OPF with `opf.ac.solver` set to `DEFAULT` with TSPOPF installed.
    *Thanks to Carlos Murillo-Sanchez.*
  - Fix bug in setup of OPF (introduced since v6.0), triggered by running a
    DC OPF with `opf.ac.solver` set one of the legacy MEX-based solvers such
    as `PDIPM`. Set `opf.ac.solver` to `PDIPM` for some of the DC OPF tests.

#### 10/10/17
  - Add `savechgtab()` function to save change tables, such as those used
    by `apply_changes()`, to a file.
  - Fix issues with PSS/E import on newer Octave versions (e.g. 4.3).
    Fixes to MATLAB incompatibilities in `regexp` behavior eliminated the
    need for Octave-specific workarounds.

#### 9/26/17
  - Minor updates to support the latest versions of MATLAB, MOSEK and YALMIP.

#### 9/15/17
  - Add another purely synthetic case from the ACTIVSg team (ASU, Cornell,
    Texas A&M, U of IL, and VCU - Synthetic grids), resulting from work
    supported by the ARPA-E GRID DATA program.
    *Thanks to Adam Birchfield and the ACTIVSg team.*
    - `case_ACTIVSg10k` (10,000-bus US WECC synthetic model)

#### 9/12/17
  - Update `@opt_model` API for method naming consistency. Summary of
    deprecated method names, with new alternatives in parenthesis:
    - `add_vars` (`add_var`)
    - `add_costs` (`add_legacy_cost`, `add_quad_cost` or `add_nln_cost`)
    - `add_constraints` (`add_lin_constraint` or `add_nln_constraint`)
    - `build_cost_params` (no longer needed)
    - `compute_cost` (`eval_legacy_cost`)
    - `get_cost_params` (`params_legacy_cost`)
    - `getv` (`params_var`)
    - `linear_constraints` (`params_lin_constraint`)

#### 9/8/17
  - Add general nonlinear cost values to OPF results in `results.nlc`.

#### 9/6/17
  - Refactor AC and DC OPF code to use the new quadratic and general
    nonlinear cost handling of `opt_model` to build and evaluate
    generator and user-defined costs and their derivatives.

#### 9/5/17
  - Add tests for OPF with high-degree polynomial (greater than quadratic)
    generator cost functions.
  - Add tests for OPF with legacy user-defined cost functions that include
    "dead zone" with quadratic "penalty".
  - Lay ground work for user-defined nonlinear OPF costs by adding
    support to `opt_model` for handling nonlinear costs with function
    handles for evaluating the cost function, gradients and Hessian.
  - Add support to `opt_model` for handling quadratic costs.
  - Deprecate the legacy generalized cost mechanism in `opt_model`
    based on `add_costs()` in favor the quadratic and general nonlinear
    mechanisms, `add_quad_cost()` and `add_nln_cost()`.

#### 8/22/17
  - Add options `'cpf.enforce_v_lims'` and `'cpf.enforce_flow_lims'` to
    enforce bus voltage magnitude and branch flow limits in the
    continuation power flow, and `'cpf.v_lims_tol'` and
    `'cpf.flow_lims_tol'` to control the respective detection tolerances.
    *Thanks to Ahmad Sadiq Abubakar and Shrirang Abhyankar.*

#### 8/18/17
  - Expand support for soft branch flow limits in `toggle_softlims` to
    include AC OPF problems as well as DC OPF.
  - Add support for direct specification of user-defined nonlinear
    constraints for AC OPF, in `mpc.user_constraints.nle` and
    `mpc.user_constraints.nli`, for equality and inequality constraints,
    respectively.

#### 8/14/17
  - Remove `nln` (nonlinear constraints) field from `opt_model` object,
    in favor of `nle` (nonlinear equalities) and `nli` (nonlinear
    inequalities).
  - Add `nle.lambda.<name>` and `nli.mu.<name>` to OPF `results` struct.
  - Add `nle` and `nli` fields to `results.mu` returned by `*opf_solver`
    functions.
  - **INCOMPATIBLE CHANGE:** Remove `nln.mu.l.<name>` and `nln.mu.u.<name>`
    fields from OPF `results` struct. Use `nle.lambda.<name>` and
    `nli.mu.<name>` fields instead for nonlinear constraint multipliers.
  - **INCOMPATIBLE CHANGE:** Modify order of default output arguments of
    `opt_model/get_idx()`.
  - **INCOMPATIBLE CHANGE:** Add `mpopt` to OPF `formulation` callback
    input args.

#### 8/4/17
  - Refactor AC OPF code to use the new nonlinear constraint handling
    of `opt_model` to build and evaluate power balance and branch flow
    constraints and their derivatives, and index shadow prices.
  - Add option for `opt_model/get_idx()` to return arbitrarily selected
    index types.

#### 7/10/17
  - Lay ground work for user-defined nonlinear OPF constraints by adding
    support to `opt_model` for handling nonlinear constraints with
    function handles for evaluating the constraint function, gradients
    and Hessian.
  - Deprecate the `add_constraints()` method of `opt_model`. Use the
    corresponding one of the following methods instead:
    `add_lin_constraint()`, `add_nln_constraint()` or `init_indexed_name()`.

#### 6/8/17
  - Move `@opt_model`, `@opf_model` to use `classdef`. Modify code to use
    OOP notation everywhere for `@opt_model`/`@opf_model` objects,
    e.g. `om.property`, `om.method()`.
    __IMPORTANT NOTE__: For Octave users, the minimum requirement for
    MATPOWER is now Octave 4 or later.

#### 5/25/17
  - Add option to call `@opt_model/compute_cost` without `idx` argument
    and have it total over all indices for a given `name`.

#### 5/24/17
  - Remove OPF result columns from `gen` matrix returned by `psse2mpc`.
  - Limit length of name of star-buses added by `psse2mpc` to 12 chars.
  - Add `save2psse` with support for exporting MATPOWER case data to
    PSS/E RAW format.

#### 5/23/17
  - Add support for `gentype` and `genfuel` fields of MATPOWER case struct
    in `extract_islands`, `ext2int`, `int2ext`, `load2disp`  and `savecase`.
  - Add support for `bus_name` field of MATPOWER case struct to
    `extract_islands`, `ext2int` and `int2ext`.
  - Add `gentype` and `genfuel` fields to three ACTIVSg cases.

#### 5/22/17
  - Add `genfuels` and `gentypes` to establish standard set of values for
    optional `mpc.genfuel` and `mpc.gentype` fields for generator fuel
    type and generator unit type, respectively.
  - Fix bug #13 where setting all buses to type `NONE` (isolated) resulted
    in a fatal error for `ext2int`, `runpf`, `runcpf` and `runopf`.
    *Thanks to SNPerkin.*

#### 5/19/17
  - Add three new purely synthetic cases from the ACTIVSg team (ASU, Cornell,
    Texas A&M, U of IL, and VCU - Synthetic grids), resulting from work
    supported by the ARPA-E GRID DATA program.
    *Thanks to Adam Birchfield and the ACTIVSg team.*
    - `case_ACTIVSg200` (200-bus Illinios synthetic model)
    - `case_ACTIVSg500` (500-bus South Carolina synthetic model)
    - `case_ACTIVSg2000` (2000-bus Texas synthetic model)

#### 5/11/17
  - Fix bug #12 where the CPF could terminate early when requesting
    trace of the full curve with P or Q limits enforced, if a limit
    becomes binding at the base case.
    *Thanks to Felix.*

#### 5/3/17
  - Fix #11 fatal error encountered when running `test_matpower` with
    SDP_PF and YALMIP installed, but no SDP solver. Now checks for
    availability of SeDuMi, SDP3 or MOSEK before attempting to run
    SDP_PF tests that require solving an SDP.

#### 4/7/17
  - Fix fatal bug in `get_losses` when computing derivatives of reactive
    branch injections and fix some related tests.
  - Fix bug in `makeJac` in which voltage was set by generator voltage
    setpoint even for PQ buses.

#### 4/7/17
  - Fix fatal bug #8 when calling `runcpf` with base and target cases with
    identical load and generation.
    *Thanks to Felix.*

#### 3/17/17
  - In the Newton power flow, for larger systems use explicit LU
    decomposition with AMD reordering and the 3 output argument form of LU
    (to select the Gilbert-Peierls algorithm), resulting in up to a 2x
    speedup in MATLAB, 1.1x in Octave.
    *Thanks to Jose Luis Marin.*
  - Add new `pf.nr.lin_solver` option to control the linear solver used
    to compute the Newton update step in the Newton-Raphson power flow.

#### 2/9/17
  - Add three new power flow algorithms for radial distribution
    systems selected via the three new options for `pf.alg`, namely
    `'PQSUM'`, `'ISUM'`, `'YSUM'`. Also includes new MATPOWER options
    `pf.radial.max_it` and `pf.radial.vcorr`. See Section 4.3 on
    "Distribution Power Flow" in the manual for details.
    *Thanks to Mirko Todorovski.*
  - Add 6 new radial distribution system cases. *Thanks to Mirko Todorovski.*
    - `case4_dist`
    - `case18`
    - `case22`
    - `case69`
    - `case85`
    - `case141`

#### 1/26/17
  - Add sources for MATPOWER User's Manual to `docs/src`.
  - Move MOST User's Manual from `docs` to `most/docs`.

#### 1/25/17
  - Move case files to new `data` directory. Requires user to update
    MATLAB or Octave path.

#### 1/24/17
  - Update documentation in `README.md`, `CHANGES.md`,
    `CONTRIBUTING.md`, `relnotes/RELEASE-NOTES-6.0.md`,
    `docs/MATPOWER-dev-guide.md` to go along with move to GitHub
    and addition of installer.

#### 1/23/17
  - Add `install_matpower()` to assist with installation by
    updating MATLAB or Octave paths or providing the commands
    required to so.

#### 1/16/17
  - Support plotting of multiple nose curves in CPF by allowing
    option `cpf.plot.bus` to take on vector values.

#### 1/14/17
  - Add line for curtailed load to `case_info()` output.

#### 1/5/17
  - Fix bug #4 where some Q limits were not being respected by CPF
    when buses were converted to PQ by initial power flow run.
    *Thanks to Shruti Rao.*

#### 1/4/17
  - When `genfuel` field is present in `mpc`, `load2disp()` now augments
    it with `dl` entries for the dispatchable loads it adds.

#### 1/3/17
  - Change default implementation of active power line flow
    constraints (`opf.flow_lim = 'P'`) to use flow directly, rather
    than square of flow, which is now a separate option, namely
    `opf.flow_lim = '2'`. *Thanks to Nico Meyer-Huebner.*

#### 12/29/16
  - Fix bug in converting older versions of MATPOWER options struct.

#### 12/21/16
  - Add Travis-CI integration. *Thanks to Richard Lincoln.*


Version 6.0 - *Dec 16, 2016*
----------------------------

#### 12/16/16
  - Released 6.0.
  - Moved development to GitHub: <https://github.com/MATPOWER/matpower>.

#### 12/9/16
  - Bumped MIPS version to 1.2.2.
  - Renamed MIPS from MATLAB Interior Point Solver to MATPOWER Interior
    Point Solver.

#### 12/6/16
  - Remove dependence of `t_mpsolve()` on presence of `have_fcn()` to
    detect PARDISO installation.
  - Remove dependence of `mpver()` on presence of `mostver()`.


Version 6.0b2 - *Nov 1, 2016*
-----------------------------

#### 11/1/16
  - Released 6.0b2.

#### 10/28/16
  - Added new general event handling mechanism to continuation power flow
    (CPF), including handling of generator real and reactive power limits.

#### 10/26/16
  - Include support for MOSEK 8, some algorithms eliminated, others changed
    algorithm code set via `mosek.lp_alg` option.

#### 10/19/16
  - Improved robustness of infeasibility detection when `pf.enforce_q_lims`
    option is used.
  - Add tests for power flow with `pf.enforce_q_lims` option.

#### 10/14/16
  - Fixed bugs in AC OPF solver code for fmincon, IPOPT, and Knitro
    involving shadow prices on variable bounds, when upper and lower bounds
    are equal (such as a voltage magnitude setpoint when `opf.use_vg` is 1).
    In IPOPT the prices are missing, and Knitro and fmincon can both return
    negative prices on the wrong constraint.

#### 10/13/16
  - Added `opf.use_vg` option to provide a convenient way to have the
    OPF respect the generator voltage setpoints specified in the gen
    matrix.

#### 10/4/16
  - Improve backward compatibility of `mpoption()`, allow it to handle
    old-style option vectors from earlier versions of MATPOWER, namely
    3.2 and 4.0, and to accept on old style options vector with new-style
    name/value pair option overrides.
  - **INCOMPATIBLE CHANGE:** Remove `cpf.user_callback_args` option and
    modify `cpf.user_callback` to allow for structs defining callback
    priority and args, in addition to just callback function name.
  - Add user options for setting tolerances for target lambda detection
    and nose point detection, `cpf.target_lam_tol` and `cpf.nose_tol`,
    respectively.

#### 9/27/16
  - **INCOMPATIBLE CHANGE:** Changed name of `cpf.error_tol` option to
    `cpf.adapt_step_tol`.

#### 9/22/16
  - Update tests for compatibility with Optimization Toolbox 7.5 (R2016b),
    which removes `active-set` algorithm for `quadprog()` and `active-set`
    and `simplex` for `linprog()`.
  - Fix fatal error when using some SDP_PF functions when CPLEX is not
    installed.

#### 9/19/16
  - Fix dimension bug in `makeBdc()` when last bus is not connected (should
    never happen if bus is properly marked as type `NONE`). *Thanks to
    Samuel Perkin.*

#### 8/15/16
  - Fix a harmless bug in `@opt_model` where variable, constraint and
    cost sets indexed by a single variable would allocate a square
    matrix of starting and ending indices, rather than a simple vector.
    *Thanks to Alberto Lamadrid for catching this.*

#### 7/18/16
  - Extract computation of tangent vector from `cpf_predictor()` into
    `cpf_tangent()`, and call sequentially as needed.

#### 7/11/16
  - Add step size to CPF callback args and results logging.
  - Add option `cpf.adapt_step_damping` to control oscillations in
    adaptive step size control for continuation power flow.

#### 7/8/16
  - Updated network reduction code to handle cases with radially connected
    external buses.

#### 7/1/16
  - Fix bug in `cpf_default_callback()` that sometimes resulted in plotting
    voltage at different buses for different parts of the curve when the
    bus was not explicitly specified with the `cpf.plot.bus` option.

#### 6/23/16
  - Bad bus numbers no longer cause a fatal error (after reporting the
    bad bus numbers) in `case_info()`.
  - Updated versions of `qcqp_opf()` and `qcqp_opf()` in `extras/misc`, from
    Cedric Josz.

#### 6/10/16
  - Fix bug in `savecase()` where single quotes were not escaped properly
    in bus names.
  - Generator capability curve parameters that define a zero-reactive
    power line no longer cause a fatal error.


Version 6.0b1 - *Jun 1, 2016*
-----------------------------

#### 6/1/16
  - Released 6.0b1.

#### 5/27/16
  - **INCOMPATIBLE CHANGE:** Removed `fairmax()` from the public interface
    by moving it inside `uopf()`, the only place it was used.

#### 5/25/16
  - Add contributed code from Camille Hamon to `extras/maxloadlim` for
    finding the loadability limits in power systems based on an optimal
    power flow using dispatchable loads.

#### 5/24/16
  - Fix bug in `toggle_dclines()` that resulted in fatal error when
    used with OPF with reactive power costs. *Thanks to Irina Boiarchuk.*
  - Add option to call `total_load()` with full case struct, instead
    of separate `bus` and `gen` matrices.

#### 5/18/16
  - Add `plot_mpc()`, contributed by Paul Cuffe, to `extras/misc`. Plots
    an electrically meaningful drawing of a MATPOWER case.
  - Add `case145.m`, IEEE 145 bus, 50 generator dynamic test case from
    http://www.ee.washington.edu/research/pstca/dyn50/pg_tcadd50.htm.

#### 5/17/16
  - Add 9 new case files, 8 cases ranging from 1888 to 6515 buses
    representing the French system, and a 13,659-bus case representing
    parts of the of the European high voltage transmission network,
    stemming from the Pan European Grid Advanced Simulation and State
    Estimation (PEGASE) project. *Thanks again to Cedric Josz and
    colleagues from the French Transmission System Operator.*
  - Add `extras/misc/qcqp_opf.m`, by Cedric Josz, et. al.

#### 5/3/16
  - Fix fatal bug in `update_mupq()` affecting cases where `QMIN` is
    greater than or equal to `QC1MIN` and `QC2MIN` (or `QMAX` is less than
    or equal to `QC1MAX` and `QC2MAX`) for all generators.
    *Thanks Jose Miguel.*

#### 3/29/16
  - Add support for `quadprog()` under GNU Octave.

#### 3/3/16
  - Copying a field containing a struct to a non-struct field with 
    `nested_struct_copy()` now overwrites rather than causing a fatal
    error.

#### 2/29/16
  - Added option to call `scale_load()` with full case struct, with
    `cost` field in `opt` to indicate when to include cost scaling.

#### 2/22/16
  - Add major new feature: *MATPOWER Optimal Scheduling Tool* (MOST).
    See `docs/MOST-manual.pdf` for details.

#### 2/18/16
  - Updated code from 9/23/16 to turn off pesky CPLEX warnings to
    included CPLEX 12.6.3.

#### 2/12/16
  - Added Release History section to Appendix of manual.

#### 1/29/16
  - Added option to call `makePTDF()`, `makeB()`, and `makeBdc()` with
    `mpc` struct instead of individual `baseMVA`, `bus`, `branch` args.
    *Suggested by Alberto Lamadrid.*

#### 1/26/16
  - Introduced work-around and warning for crashes caused by strange
    behavior from MATLAB's `ver()` function when MATPOWER (or any other
    3rd party toolbox with a `Contents.m`) is installed in a directory on
    the MATLAB path named `matlab` or `optim` (both case insensitive).

#### 1/15/16
  - Added `feval_w_path()` function for evaluating functions located at
    a specified path, outside of the MATLAB path.

#### 1/14/16
  - Added tests for `loadcase()` for m-file cases outside the MATLAB path.

#### 1/6/16
  - Added `apply_changes()` and `idx_ct()` to implement general method for
    applying modifications to an existing MATPOWER case.

#### 11/5/15
  - Use voltage dependent loads in both base and target injections
    to define continuation power flow transfer. Should fix issue with
    calculation of final loads.

#### 11/4/15
  - Fixed a bug in `psse_convert_xfmr()` where conversion of data for
    transformers with CZ=3 was done incorrectly. *Thanks to Jose Marin
    and Yujia Zhu.*

#### 10/30/15
  - Fixed a bug in `cpf_default_callback()` introduced with the
    experimental updates regarding ZIP loads on 4/14/15.

#### 10/15/15
  - Modified `t_is()` to handle matrix inputs of dimension greater
    than two.
  - Added `t_test_fcns()` to test `t_ok()` and `t_is()` and manually
    check output of failed tests.
  - Added tests to `t_dcline()` for an isolated generator bus connected
    to the rest of the system via a DC line. This works for AC and DC
    power flow and OPF, though you must set the isolated bus to
    be of type `REF`.

#### 9/23/15
  - Added code in `cplex_options()`, `insolvablepfsos()`,
    `insolvablepfsos_limitQ()` and `yalmip_options()` to turn off
    `MATLAB:lang:badlyScopedReturnValue` warning triggered by
    CPLEX when using MATLAB R2015b (8.6) and later.

#### 7/16/15
  - Added `mpopt2qpopt()` to provide common interface for creating
    options struct for `mi/qps_matpower()` from a MATPOWER options
    struct.
  - Changed default solver order for LP, QP, MILP, MIQP problems
    to move Gurobi before CPLEX and BPMPD after OT and GLPK.

#### 7/9/15
  - Modified `have_fcn()` to return 0 and warn for unrecognized
    functionality, instead of producing fatal error.

#### 6/19/15
  - Fixed a fatal bug in `psse_convert_xfmr()` affecting transformers
    with CW and/or CZ equal to 1. *Thanks to Matthias Resch.*

#### 5/18/15
  - Fixed a crash in `have_fcn()` caused by changes in OPTI Toolbox
    v2.15 (or possibly v2.12).

#### 4/24/15
  - Commented out isolated bus 10287 in `case3375wp.m`.

#### 4/16/15
  - Added some caching to `mpoption()` and made minor changes to
    `nested_struct_copy()` to greatly decrease the overhead added by
    `mpoption()` when running many small problems.

#### 4/14/15
  - Added experimental code to lay foundation for handling ZIP load
    model in power flow (Newton, fast-decoupled only), continuation
    power flow, and optimal power flow (MIPS, fmincon, Knitro, IPOPT
    solvers only). Currently, ZIP loads can only be specified on a
    system-wide basis using the experimental options
    `exp.sys_wide_zip_loads.pw` and `exp.sys_wide_zip_loads.qw`.
    Tests in `t/t_vdep_load()`.
  - Added `bus` and `area` as possible values for `load_zone` argument
    to `total_load()`, which is now used to compute voltage dependent
    load values.

#### 3/27/15
  - Fixed issue where default value of `feastol` option was not being
    set correctly in `mips()` when called directly (or via `qps_mips()`)
    with `feastol = 0`. By default, MATPOWER option `mips.feastol` is
    set to zero but is normally replaced in `mipsopf_solver()` or
    `dcopf_solver()` with value of `opf.violation` option before
    calling `mips()`, thereby masking the problem.
  - In `miqps_glpk()` variables of type `B` (binary) are converted to
    `I` (integer) with appropriate bounds, since some versions of
    GLPK do not natively handle type `B`.

#### 3/24/15
  - Added code to DC OPF to return `success` = 0 for cases where the
    matrix is singular (e.g. islanded system without slack).
  - Fixed problem in `have_fcn()` where SeDuMi was turning off and
    leaving off all warnings.


Version 5.1 - *Mar 20, 2015*
----------------------------

#### 3/20/15
  - Released version 5.1.

#### 3/19/15
  - Added support for using PARDISO as linear solver for computing
    interior-point update steps in MIPS (v1.2), via new `mplinsolver()`
    function and `mips.linsolver` option.

#### 3/18/15
  - Added four new case files, ranging from 89 up to 9421 buses,
    representing parts of the European high voltage transmission
    network, stemming from the Pan European Grid Advanced Simulation
    and State Estimation (PEGASE) project. *Thanks to Cedric Josz and
    colleagues from the French Transmission System Operator.*
  - Added network reduction toolbox to `extras/reduction` directory
    for creating smaller approximate network equivalents from a larger
    original case file. *Thanks to Yujia Zhu and Daniel Tylavsky.*

#### 3/4/15
  - Added variable type as an attribute to `@opt_model`, so you can
    now specify variables as `C`, `I`, or `B` (continuous, integer,
    or binary) when adding variables with `add_vars()` and `getv()` can
    optionally return a variable-type string suitable for passing
    to `miqps_matpower()` and friends.

#### 2/26/15
  - Minor speed improvements in various `@opt_model` functions from
    bypassing calls to `substruct()`.

#### 2/25/15
  - Switch to more permissive 3-clause BSD license from GPL 3.0.

#### 2/24/15
  - Added function `mpoption_info_intlinprog()`, tag `intlinprog` to
    `have_fcn()` and optional `intlinprog` field to MATPOWER options.
  - Added `miqps_matpower()`, a wrapper function for various solvers
    of mixed-integer linear and quadratic programs. Functionality
    implemented by `miqps_cplex()`, `miqps_glpk()`, `miqps_gurobi()`,
    `miqps_mosek()` and `miqps_ot()`. Added corresponding
    `t/t_miqps_matpower()` to test suite.

#### 2/18/15
  - Changed generator and dispatchable load sections in `printpf()`
    output to include off-line units, after all, it already has a
    Status column.

#### 2/13/15
  - Added explicit colors to plots in `cpf_default_callback()` so things
    look right in newer MATLAB versions (R2014b and later).

#### 2/6/15
  - Modified `nested_struct_copy()` to eliminate `cellfun()` call that
    was not supported in MATLAB 7.0.x. Noted that `runcpf()` also
    requires MATLAB 7.1 or later due to a call to `cellfun()`.
    Included code to skip certain tests that require that `cellfun()`
    functionality when running under MATLAB 7.0.x.

#### 2/5/15
  - Replaced `regexp(... 'split')` construct in `mpoption()` with
    `regexp(... 'tokens')` since it was causing fatal errors on
    MATLAB versions < 7.3, which did not have that feature.
  - Fixed fatal error in when using fast-decoupled power flow
    on MATLAB versions < 7.3, caused by use of newer
    `lu(... 'vector')` construct.

#### 2/3/15
  - Added check to `have_fcn()` for installation of `ipopt_auxdata.m`
    when Ipopt >= 3.11.x is detected, to warn about incomplete
    installation and avoid a fatal error.

#### 1/27/15
  - Added online function reference, produced by m2html. *Thanks to
    Shrirang Abhyankar.*

#### 1/26/15
  - Fixed bug in `hasPQcap()` that resulted in ignoring generator
    capability curves if `Q1MAX < Q2MAX` or `Q1MIN > Q2MIN` (i.e. when
    the sloped portions cut off the left corners of the box
    constraints, rather than the right corners).
    *Thanks to Irina Boiarchuk.*

#### 1/23/15
  - Added `mosek_symbcon()` to define symbolic constants for setting
    MOSEK options. Updated MOSEK solver option values in help text
    for `mpoption()` to correspond to MOSEK v7.

#### 1/22/15
  - Added ability to toggle the availability of optional functionality
    using `have_fcn()`.

#### 1/21/15
  - Fixed minor bug with `poly2pwl()`, affecting units with `Pmax <= 0`.
  - Major update to `have_fcn()`, which now determines and caches
    version numbers and release dates for optional packages, used
    by `mpver()` and others.
  - Fixed error in `qps_mosek()` in printout of selected optimizer
    when using MOSEK 7.

#### 1/16/15
  - Cleaned up and improved consistency of output in `printpf()` for
    generation and dispatchable load constraints.
  - Modified `runcpf()` to gracefully handle the case when the base
    and target cases are identical (as opposed to getting lost in
    an infinite loop).

#### 1/15/15
  - Fixed bug in handling of interface flow limits, where multipliers
    on binding interface flow limits were off by a factor of the p.u.
    MVA base.
  - Fixed sign error on multipliers on lower bound on constraints
    in `qps_clp()` and `qps_glpk()`. Modified a test in `t_qps_matpower()`
    to check for this.

#### 1/14/15
  - Added support for LP/QP solver CLP (COIN_OR Linear Programming).
    Use `opf.dc.solver` option `CLP` or `qps_clp()`.
  - Added note to README and home page about OPTI Toolbox by
    Jonathan Currie being an easy way to install some good solvers
    for Windows users.

#### 1/13/15
  - Updated `t_opf_dc_ot()` to remove the skipping of the checking of
    price results for the dual-simplex algorithm for all versions of
    MATLAB except R2014b, the first version that included the
    dual-simplex algorthm. For some reason, in this version it did
    not return any Lagrange multipliers!?!.
  - Increment MATPOWER options version number to 5 (forgot to do it
    previously for 3, 4, and 5) and included code to update older
    options structs to current version.

#### 1/7/15
  - Improved detection of GLPK version in `mpver()` and GLPK
    availability in `have_fcn()`, now compatible with GLPK installed
    by OPTI Toolbox (http://www.i2c2.aut.ac.nz/Wiki/OPTI/).

#### 12/22/14
  - Fixed fatal bug in `toggle_dcline()` when pretty-printing results.
    *Thanks to Deep Kiran for reporting.*

#### 12/18/14
  - Fixed fatal bug in `case_info()` for islands with no generation.


Version 5.0 - *Dec 17, 2014*
----------------------------

#### 12/17/14
  - Released version 5.0.

#### 12/16/14
  - Added unsupported functions `check_feasibility.m`, `checklimits.m`,
    `loss2bus.m`,  `make_opf_feasible.m`, and `makeBloss.m` to `extras/misc`.
  - Added section 9 "Miscellaneous MATPOWER Functions" to User's Manual.

#### 12/15/14
  - Added option for `case_info()` to print to a file. Added `case_info()`
    tests to `t_island()`.

#### 12/12/14
  - Added code in `psse_read()` to work around a bug in MATLAB 7.3.
  - Added private `catchme` tag (for internal use only) to `have_fcn()`
    to detect older versions of MATLAB and Octave that do not support
    the `catch me` syntax in try/catch constructs.
  - Added private `regexp_split` tag (for internal use only) to
    `have_fcn()` to detect older versions of Octave that do not support
    the `split` argument to `regexp()`.
  - Updated to support Ipopt 3.11.x and later, which removed support
    for `options.auxdata` from the MEX file. Added private
    `ipopt_auxdata` tag (for internal use only) to `have_fcn()` to
    detect version 3.11.x or later.

#### 12/4/14
  - Added new option `opf.init_from_mpc` to force some solvers to
    use the starting point supplied in the MATPOWER case to
    initialize the optimization variables for the OPF, instead
    of creating its own starting point. Currently only implemented
    for Ipopt, Knitro and MIPS.
  - **INCOMPATIBLE CHANGE:** Renamed `cdf2matp()` to `cdf2mpc()` and modified
    the interface to be consistent with `psse2mpc()`.

#### 12/2/14
  - Added new option `out.suppress_detail` to quickly suppress all
    pretty-printed output except the system summary. By default,
    detailed output is automatically suppressed for systems larger
    than 500 buses.
  - DC OPF formulation in `opf_setup()` now uses a single set of
    doubly-bounded constraints for flow limits, instead of two
    sets of upper bounded constraints.
  - Updated to MIPS 1.1, which includes additional user-settable
    options: `xi`, `sigma`, `z0`, `alpha_min`, `rho_min`, `rho_max`,
    `mu_threshold` and `max_stepsize`.
  - **INCOMPATIBLE CHANGE:** The name of the option to `mips()` to specify
    the maximum number of step-size reductions when `step_control` is on
    was changed from `max_red` to `sc.red_it` for consistency with
    other MATPOWER options.

#### 11/18/14
  - Updated `case300.m` with new conversion from original CDF file.
    No longer uses 9900 MVA as proxy for unlimited line capacity.

#### 11/13/14
  - Added capability for `set_reorder()` to automatically pad matrix
    or cell array before assignment of matrix or cell array with
    larger size in some dimensions. Added tests to `t_ext2int2ext()`.

#### 11/12/14
  - Loads at isolated buses are no longer included in results from
    `total_load()`.
  - Fixed loads and shunts at isolated buses are no longer included
    in system and area summaries in `printpf()`. Line charging
    injections set to zero for DC power flow results.
  - Added tests for `printpf()` in `t_printpf.m`.
  - Changes to the internally indexed version of `gencost` now
    get properly copied back to externally indexed version by
    `int2ext()`.
  - Added tests in `t_ext2int2ext()` to confirm that `e2i_data/field()`
    and `i2e_data/field()` work for cell array as well as numerical
    array fields. Updated help text to reflect this feature.

#### 11/11/14
  - Added `get_losses()` function to compute losses, line charging
    reactive injections and their derivatives, as functions of
    bus voltages. Corresponding tests included in `t_get_losses()`,
    including example of loss sensitivity factors.
  - Fixed bug in `runpf.m` that caused a crash for cases with
   ` pf.enforce_q_lims` turned on and exactly two Q limit violations,
    one Qmax and one Qmin. *Thanks to Jose Luis Marin.*

#### 11/4/14
  - Modified behavior so setting `out.all` option to 0 is now ignored
    for pretty-printed output to files specified as `FNAME` argument
    to `runpf()`, `runopf()`, etc.

#### 10/22/14
  - Removed `idx_area()` and all code references to `areas` field
    in MATPOWER cases except those needed to support reading
    v1 case files and those required for backward compatibility of
    APIs. Removed unused (and formerly deprecated) `areas` field from
    version 2 case files that still included it.

#### 10/21/14
  - Fixed a bug in `savecase()` where a `gencost` matrix with extra
    columns of zeros resulted in a corrupted MATPOWER case file.
  - Modified `total_load()` to return actual rather than nominal
    value for dispatchable loads by default, unless using the
    old-style string input for the 4th input arg. See
    `help total_load` for details.

#### 10/16/14
  - Reactive power output of multiple generators at a PQ bus
    no longer get re-allocated when running a power flow.

#### 10/11/14
  - Added `fmincon_ip`, `linprog_ds` and `optimoptions` to
    `have_fcn()` to test for availability of fmincon's interior
    point method (Optimization Toolbox 4.x +), linprog's
    dual simplex method (Optimization Toolbox 7.1 +), and
    optimoptions function for setting Optimization Toolbox
    options (Optimization Toolbox 6.3 +), respectively.
  - Added handling of NaN values to `t_is()`.
  - **INCOMPATIBLE CHANGE:** Removed use of `ot_opts` field and
    replaced with `linprog_opts` and `quadprog_opts` fields
    in the `OPT` argument to `qps_matpower()` and `qps_ot()`.
  - Added optional `linprog` and `quadprog` fields to MATPOWER
    options struct, to allow setting of their options directly
    from `mpoption()`. Incremented MATPOWER options struct version
    to 2.
  - Added tests for multiple `linprog()` algorithms to `t_opf_dc_ot()`.

#### 10/10/14
  - Added `sqp` algorithm option for `fmincon.alg`, not suitable
    for large problems.

#### 10/7/14
  - Added `toggle_softlims()`, an extension to implement DC OPF
    branch flow soft limits. This should be useful in identifying
    the cause of infeasibility in some infeasible DC OPF problems.
    See `help toggle_softlims` for details.

#### 9/24/14
  - Fixed issue with failed tests in `t_psse()` on Windoze.

#### 8/11/14
  - Modified `savecase()` to save fields `bus_name` and `genfuel`
    (both cell arrays of strings) if present. This is an
    experimental feature and the design may change in future
    versions.
  - Optionally return success flag from `t_ok()` and `t_is()`.

#### 8/4/14
  - Additional improvements to correctly handling PSS/E RAW
    files for revisions 24-28.

#### 7/30/14
  - Fixed a bug in `psse_convert()` that resulted in incorrect
    bus voltage angles when importing from v29 or v30 files.
  - Enhanced PSS/E import code to handle additional versions
    prior to v29, improve reporting of version being used
    for parsing. Default version when not explicit is now v23.

#### 7/24/14
  - Added line to `cplex_options()` to prevent default creation
    of `clone1.log` files when using parallel CPLEX routines
    (only effective for CPLEX 12.6.0.1 or later).
  - Fixed fatal error when `uopf()` shuts down all gens
    attempting to satisfy Pmin limits.

#### 7/23/14
  - Fixed error in User's Manual description of LODF and
    further clarified text.


Version 5.0b1 - *Jul 1, 2014*
-----------------------------

#### 7/1/14
  - Released version 5.0b1.

#### 6/27/14
  - Added support for LP solver GLPK which is built-in to
    Octave. Use `opf.dc.solver` option `GLPK` or `qps_glpk()`.
  - **INCOMPATIBLE CHANGE:** Removed optional `max_it` field from
    opt argument to `qps_matpower()` and `qps_*()` family of
    functions (except `qps_mips()`).

#### 6/13/14
  - Added ability for `loadcase()` to load MATPOWER case M-files
    that are not in the MATLAB path by specifying an explicit
    path in the input filename.

#### 6/9/14
  - Fixed memory issue resulting from nested om fields when
    repeatedly running an OPF using the results of a previous
    OPF as input. *Thanks to Carlos Murillo-Sanchez.*

#### 5/23/14
  - Added `case5.m`, a 5-bus, 5-generator example case from Rui Bo.
  - Removed `extras/psse2matpower`.
  - Fixes to `printpf()` to suppress printing of dispatchable load
    constraint section when there are no dispatchable loads and
    to never print line constraints for unconstrained lines.
  - Fixed crash when using Knitro to solve cases with all
    lines unconstrained.
  - Switched to `[L,U,p,q] = lu(B,'vector')` form for factorization
    in fast-decoupled power flow in `fdpf()`, resulting in ~5x
    speedup on some large systems. *Thanks to Carlos Murillo-Sanchez.*

#### 5/6/14
  - Modified `savecase()` to automatically add a comment line
    with the function name, and make sure the name is converted
    to a legal function name before saving.
  - Further additions to PSS/E import code to handle repeated
    delimiters, accumulation of warnings to save in comments
    of converted file, improved consistency of verbose option,
    and some automated tests.

#### 4/28/14
  - Major revision of PSS/E import code to further improve
    robustness and include Octave support. Updated `psse2mpc()` to
    include direct saving of result to MATPOWER case file.
    - Added functions:
      - `psse_parse()`
      - `psse_parse_line()`
      - `psse_parse_section()`
    - Removed functions:
      - `psse_count_lines()`
      - `psse_extract_data()`
      - `psse_read_section()`

#### 4/11/14
  - Significant updates to PSS/E import code to improve robustness
    and (mostly) handle versions 29 to 33. Renamed the 4 functions
    ending in `_33` by removing the `_33`.

#### 4/7/14
  - Added experimental feature, via function `psse2mpc()`, to
    import PSS/E RAW data (version 33) into a MATPOWER case file.
    Supporting functions include:
    - `psse_convert_33()`
    - `psse_convert_hvdc_33()`
    - `psse_convert_xfmr_33()`
    - `psse_count_lines()`
    - `psse_extract_data()`
    - `psse_read_33()`
    - `psse_read_section()`

#### 3/28/14
  - Enhanced `extract_islands()` to handle DC lines, custom fields
    and extraction of multiple islands into a single case struct.

#### 3/11/14
  - Fixed bugs in `runpf()` related to enforcing generator reactive
    power limits when all generators violate limits or when
    the slack bus is converted to PQ.

#### 3/10/14
  - Fixed bug in `qps_gurobi()` where return status values for
    `NUMERIC` and `SUBOPTIMAL` were swapped. Added `INPROGRESS` status.
    *Thanks to Alberto Lamadrid for catching this.*

#### 3/2/14
  - Fixed bug in `case_info()` that incorrectly included dispatchable
    loads in amount reported for min/max capacity for Generation.

#### 2/28/14
  - Updated `toggle_dcline.m` to work correctly for OPF cases with
    user supplied constraints and costs.

#### 2/27/14
  - Fixed bug in `savecase()` where the function name mistakenly
    included the path when the `FNAME` input included a path.
    
#### 2/6/14
  - Small tweak in `connected_components()` results in ~30x
    speedup for 62k bus network. (Reminder: Never access
    indexed rows of a large sparse matrix, always transpose
    and index the columns).

#### 2/4/14
  - Added check in `qps_cplex()` for undocumented `exitstatus`
    values returned by `cplexlp()` or `cplexqp()`.

#### 1/17/14
  - Added Dan Molzahn's SDP_PF package, a set of applications of
    a semidefinite programming relaxation of the power flow
    equations, to the directory `extras/sdp_pf`.

#### 1/16/14
  - Added `status` option for 2nd argument to `toggle_reserves()`
    `toggle_dcline()` and `toggle_iflims()` as a convenient way to
    check the enabled/disabled status of these sets of callback
    functions.
  - Removed `extras/cpf` since CPF is now part of the core.
  - Added support for Dan Molzahn's SDP_PF package (coming soon).

#### 1/15/14
  - Modified handling of options for optional packages, added:
    - `mpoption_info_cplex()`
    - `mpoption_info_fmincon()`
    - `mpoption_info_gurobi()`
    - `mpoption_info_ipopt()`
    - `mpoption_info_knitro()`
    - `mpoption_info_mosek()`

#### 1/8/14
  - Updates to `qps_cplex()` and `cplex_options()` to fix verbose
    display issues with CPLEX 12.6.

#### 1/3/14
  - Added persistent variable to improve performance of `have_fcn()`.
  - Added support for Knitro v9.0.0, including new `knitromatlab`
    and `ktrlink` options to `have_fcn()`, to determine which Knitro
    interfaces are available.

#### 12/10/13
  - New MATPOWER options implementation based on options struct
    instead of options vector.
    **INCOMPATIBLE CHANGE:** In `results` struct returned by an OPF, the
    value of `results.raw.output.alg` is now a string, not an old-style
    numeric alg code.
  - Gurobi now takes precendence over MOSEK when default solver
    is selected for DC OPFs or `qps_matpower()`.

#### 12/4/13
  - Corrected error in Figure 6-5 "Total Cost Function for Negative
    Injection" in Dispatchable Loads section of User's Manual
    (slopes were labeled incorrectly).

#### 11/15/13
  - Added `case_info()`.
  - Modified `connected_components()` to sort returned groups by
    decreasing cardinality.

#### 11/5/13
  - Fixed a bug in MIPS where a near-singular matrix could produce
    an extremely large Newton step, resulting in incorrectly satisfying
    the relative feasibility criterion for successful termination.
  - Improved the starting point created for Ipopt, Knitro and MIPS by
    `dcopf_solver()`, `ipoptopf_solver()`, `ktropf_solver()` and `mipsopf_solver()`
    for variables that are only bounded on one side.

#### 10/11/13
  - Removed support for MATLAB 6.x. Removed `anon_fcns` option from
    `have_fcn()`. Files removed:
     - `fmincopf6_solver.m`
     - `mips6.m`
     - `mips6opf_solver.m`
     - `qps_mips6.m`
     - `t_mips6.m`
  - Removed support for `constr` and successive LP-based OPF solvers.
    Removed `constr`, `lp`, `qp` options from `have_fcn()`. Files removed:
     - `copf_solver.m`
     - `fun_copf.m`
     - `grad_copf.m`
     - `LPconstr.m`
     - `LPeqslvr.m`
     - `lpopf_solver.m`
     - `LPrelax.m`
     - `LPsetup.m`
     - `mp_lp.m`
     - `mp_qp.m`
     - `t/t_opf_constr.m`
     - `t/t_opf_lp_den.m`
     - `t/t_opf_lp_spf.m`
     - `t/t_opf_lp_spr.m`

#### 9/20/13
  - Added continuation power flow, `runcpf()`, with tangent
    predictor and Newton method corrector, *based on code
    contributed by Shrirang Abhyankar and Alex Flueck.*

#### 9/5/13
  - Fix in `smartmkt()` to avoid crash following non-convergent
    `uopf` when `mkt.lim.P.max_offer` is not defined.

#### 9/3/13
  - Fixed bug (typo) in `auction()` that could affect cases with
    a lower limit on the cleared bid price.

#### 7/30/13
  - Extended `modcost()` to optionally accept a vector of shift/scale
    factors, instead of just a scalar.

#### 6/7/13
  - Made non-convergent results more obvious by not printing
    the standard output tables (can be overridden with new
    `OUT_FORCE` option) and making the "did not converge"
    more prominent.

#### 5/1/13
  - Changed behavior of branch angle difference limits so that
    0 is interpreted as unbounded only if both `ANGMIN` and `ANGMAX`
    are zero. Added note about this to docs in various places.

#### 4/26/13
  - DC OPF now correctly sets voltage magnitudes to 1 p.u.
    in results.

#### 3/29/13
  - Performance optimizations in `@opt_model` for cases with
    large numbers of variable sets and linear constraints
    or costs specified as entire rows (all columns as
    opposed to specific var sets).
    
#### 3/13/13
  - Added to `scale_load()` the option to scale the `gencost`
    (specifically the quantity axis of the marginal cost function)
    corresponding to any modified dispatchable load. Simply add
    `gencost` as additional input and output args.
  - Empty `got` and `expected` arguments to `t_is()` now
    count as a passing test instead of an error, as long as
    the dimensions match.

#### 2/12/13
  - Fixed bug causing value of `opt.verbose` to be ignored in
    `qps_gurobi()`.

#### 12/14/12
  - Removed code in `fmincopf()` that attempts to find an interior
    starting point when using the interior point solver. It did
    not seem to help and caused errors for certain cases with
    DC lines (and probably other extensions).

#### 10/1/12
  - Updates to `have_fcn()` and `mpver()` to better handle case where
    Optimization Toolbox is installed, but with no valid license.

#### 8/30/12
  - Major speed-up in `@opt_model/linear_constraints()` by building
    transpose of `A` (assigning to full columns) then transposing
    back as opposed to building `A` directly (assigning full rows).

#### 8/1/12
  - Added function `margcost()` for computing the marginal cost of
    generation.

#### 7/20/12
  - Added utility function `@opt_model/describe_idx()` to identify
    variable, constraint or cost row indices to aid in debugging.

#### 7/18/12
  - Made `N` optional field (default is identity matrix) in
    `@opt_model/add_costs()`.
  - Added missing optional 2nd arg to `@opt_model/build_cost_params()`.

#### 6/26/12
  - Fixed a bug in the new `@opt_model/add_vars()` when adding a var
    set of dimension zero.

#### 6/18/12
  - Updated Gurobi interface for compatibility with native MATLAB
    support introduced in Gurobi 5.
    **INCOMPATIBLE CHANGE:** No longer works with older Gurobi 4.x/
    `gurobi_mex()` interface.
  
#### 5/3/12
  - Reimplementated `@opf_model` class as sub-class of the new
    `@opt_model` class, which supports indexed named sets of
    variables, constraints and costs.

#### 5/2/12
  - In `opf_setup()`, take magnitude of initial voltages at generator
    buses from `bus` matrix (`VM`), not `gen` matrix (`VG`).

#### 4/30/12
  - Fixed a bug in `int2ext()` where converting a case to internal
    ordering before calling `runpf()` or `runopf()` could result in
    a fatal error due to mismatched number of columns in internal
    and external versions of data matrices. *Thanks to Nasiruzzaman
    and Shiyang Li for reporting and detailing the issue.*
  - Fixed fatal bug in MIPS for unconstrained, scalar problems.
    *Thanks to Han Na Gwon. Bumped MIPS version to 1.0.1.*

#### 3/14/12
  - Fixed a bug in `runpf()` where it was using the wrong initial
    voltage magnitude for generator buses marked as PQ buses. Power
    flow of solved case was not converging in zero iterations as
    expected.

#### 2/29/12
  - Added a tolerance for detecting violated Q limits in `runpf()`
    when `ENFORCE_Q_LIMS` is true. *Suggested by Hongxing Ye.*

#### 1/31/12
  - Added utilities to help in working with networks with islands,
    `find_islands()`, `extract_islands()` and `connected_components()`
    and corresponding test file `t/t_islands()`.

#### 1/24/12
  - Added option to `makeJac()` to return full Jacobian instead of
    reduced version used in Newton power flow updates.

#### 1/16/12
  - Added new function `gurobiver()` for retreiving/printing Gurobi
    and Gurobi_MEX version information, since it is used multiple
    places.

#### 1/10/12
  - Moved the building of Ybus outside the reactive limit
    enforcement loop in `runpf()`. *Suggested by Shiyang Li.*

#### 1/4/12
  - Running a power flow for a case with DC lines but no `gencost`
    no longer causes an error.

#### 12/14/11
  - Modified `t/t_opf_fmincon.m` to use active-set method for testing
    `fmincopf` for MATLAB versions 7.6-7.9, since fmincon's interior
    point solver (now default) was not accurate enough in these
    versions.


Version 4.1 - *Dec 14, 2011*
---------------------------

#### 12/14/11
  - Released version 4.1.
  - Fixed bug in check for `ENFORCE_Q_LIMS` with multiple slacks
    in `runpf()`.
  - Moved printing of power flow solver name into `runpf()` so it
    doesn't get repeated when `ENFORCE_Q_LIMS` is on.

#### 12/9/11
  - Fixed problem with `qps_cplex()` when H matrix is not
    perfectly numerically symmetric.

#### 12/8/11
  - Added basic DC line modeling capability. See help for
    `toggle_dcline()` for details.
    
#### 12/1/11
  - Removed deprecated functions in `@opf_model`, `gen_lin_N()`,
    `get_nln_N()`, `get_var_N()`, use `getN()` instead.
  - Removed all references to deprecated option `OUT_RAW`.

#### 11/15/11
  - Fixed a crashing bug in computation of quadratic user-defined
    costs. *Thanks to Stefanos Delikaraoglou.*

#### 11/11/11
  - Changed default DC OPF/LP/QP solver precedence in 
    `dcopf_solver()` and `qps_matpower()` to the following:
    CPLEX, MOSEK, Gurobi, BPMPD, Opt Tbx, MIPS.
  - Minor enhancements to `cdf2matp()`, including saving of bus
    names. *Thanks to Alvaro Jaramillo Duque.*

#### 11/09/11
  - Refactored `ext2int()` and `int2ext()` into additional functions
    `e2i_field()`, `e2i_data()`, `i2e_field()` and `i2e_data()` to
    clean things up and prepare for the ability to automatically
    re-order data in cell arrays.

#### 10/31/11
  - Added three case files, all more recent variations of the
    Polish system: `case3012wp.m`, `case3120sp.m` and `case3375wp.m`.
    *Thanks to Roman Korab <roman.korab@polsl.pl>.*

#### 9/28/11
  - Increased threshold value used to filter constraint shadow
    prices in `printpf()`.
  - In `savecase()` increased precision of values saved in M-file
    case files

#### 9/16/11
  - Fixed that `qps_cplex()` would not print progress even with
    `VERBOSE > 0` with CPLEX 12.3.

#### 7/26/11
  - Fix in `qps_cplex()` for changed sign on multipliers with
    CPLEX 12.3 vs CPLEX 12.2.

#### 7/25/11
  - Fixed bug in `compare_case()` that would cause some column
    names for the `branch` matrix to be incorrectly reported.

#### 7/7/11
  - In `scale_load()`, when no load_zone is specified, it no longer
    incorrectly misses Q-only load buses.

#### 7/5/11
  - Added support for the Gurobi optimizer for large-scale linear
    and quadratic programming. To use Gurobi for the DC OPF, set
    `OPF_ALG_DC` = 700. Gurobi's various solvers can be selected via
    MATPOWER's `GRB_METHOD` option. Requires the Gurobi libraries
    available from http://www.gurobi.com/ and the Gurobi MEX
    interface available from http://www.convexoptimization.com/
    wikimization/index.php/Gurobi_mex.
  - Added function `qps_gurobi()` for solving QP and LP problems using
    the common QP solver interface used in MATPOWER. The `qps_matpower()`
    function also includes the option to use Gurobi.
  - Changed order of precendence of installed DC OPF solvers
    (i.e. LP/QP solvers) to Gurobi, MOSEK, CPLEX, BPMPD, then MIPS.

#### 6/29/11
  - Updated `t_is()` to properly print when result includes NaNs.

#### 6/17/11
  - Changed `FMC_ALG` option default to 4, fmincon defaults to
    using an interior-point method with user-supplied Hessians
    instead of an active set method.
  - Added support for the KNITRO optimization solver for large
    scale non-linear problems. Use `OPF_ALG = 600` for AC OPF.
    Requires the Optimization Toolbox from The MathWorks and
    the KNITRO libraries, available from http://www.ziena.com/.

#### 6/16/11
  - Complete rewrite of `update_mupq.m`. Should fix problems
    caused by non-zero multipliers on non-binding generator
    limits.

#### 5/17/11
  - Updated `runpf()` and `pfsoln()` to properly handle slack for
    power flow cases with islands and multiple reference buses.
    (Note: This does not include the case where `ENFORCE_Q_LIMS`
    results in temporarily converting a reference bus to a
    PQ bus and automatically finding a new suitable slack bus.)

#### 5/16/11
  - Pretty printed output from `printpf()` now includes a `*` after
    the voltage angle in the bus section for reference buses.

#### 3/31/11
  - Default value for ramp rates for dispatchable loads now
    set to `Inf` in `load2disp.m`.

#### 3/18/11
  - Fixed bug in `toggle_reserves.m` that computed the prices in
    `results.reserves.prc` incorrectly.


Version 4.0 - *Feb 7, 2011*
---------------------------

#### 2/16/11
  - Oops! Neglected to include Carlos as co-author on User's Manual.
    Updated 4.0 distribution and web-site with corrected `manual.pdf.`
    *Sorry about that Carlos!*

#### 2/7/11
  - Released version 4.0.

#### 1/18/11
  - Added `quadprog_ls` option to `have_fcn()` to check for availability
    of version of `quadprog()` with large scale solver.
  - Modified `qps_ot()` to set default options based on capabilities
    of version of Optimization Toolbox.

#### 12/16/10
  - Fixed bug in `qps_cplex()` where an infeasible problem resulted
    in a fatal error.
  - Fixed bug in `qps_mosek()` where exit flag was indicating success
    for infeasible solutions.
  - Fixed bug in `dcopf_solver()` where an infeasible problem found
    by CPLEX would result in a fatal error.


Version 4.0b5 - *Dec 13, 2010*
------------------------------

#### 12/13/10
  - Released version 4.0b5.

#### 12/2/10
  - Added to `opf_setup()` a better check on when specified generalized
    cost parameters are compatible with DC OPF.
  - Improved output of `t_is()`. Includes only elements violating
    tolerance.

#### 11/30/10
  - Fixed bug in `opf_execute()` related to automatic conversion of
    single-block piecewise linear costs to linear polynomial costs.
    Affected indexing of the `Va`, `Vm`, `Pg` and `Qg` portions of
    `results.x` and `raw.xr`.

#### 11/24/10
  - Added support for the MOSEK optimizer for large-scale linear and
    quadratic programming. To use MOSEK for the DC OPF, set
    `OPF_ALG_DC = 600`. Specific LP algorithms can be selected by
    the `MOSEK_LP_ALG` option. Requires the MATLAB interface for MOSEK,
    available from http://www.mosek.com/.
  - Added function `qps_mosek()` for solving QP and LP problems using
    the common QP solver interface used in MATPOWER. The `qps_matpower()`
    function also includes the option to use MOSEK.

#### 11/16/10
  - Fixed bug in `opf_setup()` where indexing data for branch angle
    difference limits was not being saved in the case of DC OPF.
  - Added support for the IBM ILOG CPLEX optimizer for
    large scale linear and quadratic programming. To use CPLEX
    for the DC OPF, set `OPF_ALG_DC = 500` and choose the specific
    CPLEX solver using options `CPLEX_LPMETHOD` and `CPLEX_QPMETHOD`.
    Requires the MATLAB interface for CPLEX, available from
    http://www.ibm.com/software/integration/optimization/cplex-optimizer/.
  - Added function `qps_cplex()` for using CPLEX to solve QP and LP
    problems using the common QP solver interface used in MATPOWER. The
    `qps_matpower()` function also includes the option to use CPLEX.
  
#### 11/9/10
  - Fixed an indexing bug in `dcopf_solver()` affecting cases with a mix
    of piecewise linear and polynomial costs (unless the polynomial
    costs came after all of the pwl costs).

#### 10/12/10
  - Performance optimization in `opf_consfcn()`. Assign sparse cols
    then transpose instead of assigning sparse rows. Results in >2x
    speed up for this function on `case2935`, ~10x on `case42k`.

#### 7/20/10
  - Made some updates to `extras/psse2matpower`. Added/fixed some comments,
    text output, switched to `Parse::Text::parse_line` for bus data to fix
    problem caused by certain characters (e.g. `/` `,`) in bus names. Fixed
    error in switched shunt data (using data from wrong column). Modified
    to no longer comment out isolated buses, since it doesn't remove
    corresponding gens/branches.

#### 6/29/10
  - Fixed bug in `uopf()`. Was not calling `printpf()` when called with no
    output arguments. *Thanks to V. Ravikumar Pandi.*

#### 6/25/10
  - Added `makeJac()`, a utility function to form the power flow Jacobian.
  - Modified `makeYbus()` to allow for single MATPOWER case struct as input.
  - Added `load2disp()` to convert fixed loads to dispatchable loads.

#### 6/1/10
  - Added `modcost()` and tests for `modcost()` and `totcost()`.


Version 4.0b4 - *May 21, 2010*
------------------------------

#### 5/21/10
  - Released version 4.0b4.

#### 5/18/10
  - Added support for the IPOPT interior point optimizer for
    large scale non-linear optimization. Use `OPF_ALG = 580`
    and `OPF_ALG_DC = 400` for AC and DC OPF, respectively. Requires
    the MATLAB MEX interface for IPOPT, available from
    http://www.coin-or.org/projects/Ipopt.xml.

#### 5/13/10
  - Modified input args for Hessian evaluation function for MIPS.
    Requires `cost_mult` as 3rd argument.
  - Added check for invalid `gencost` `MODEL` in `opf_setup()`.

#### 5/5/10
  - Added `RETURN_RAW_DER` option to control when OPF solver returns
    constraint, Jacobian and objective function gradient and Hessian
    information.

#### 5/4/10
  - Refactored portions of `opf()` into `opf_setup()` and `opf_execute()`.


Version 4.0b3 - *Apr 19, 2010*
------------------------------

#### 4/19/10
  - Released version 4.0b3.
  - Changed licensing to GNU General Public license. See `LICENSE` and
    `COPYING` files for details.
  - Added in `extras` sub-directory `psse2matpower` Perl script for
    converting PSS/E data files to MATPOWER case files.

#### 4/6/10
  - Added `anon_fcns` to `have_fcn()` to check for anonymous function
    capability to avoid direct MATLAB version checks in code.
  - GNU Octave compatibility!  (tested with Octave 3.2.3)
    Added `octave` to `have_fcn()` to check when code is running under
    Octave instead of MATLAB.

#### 3/23/09
  - Fixed bug in smart market code that caused it to die for cases with
    non-consecutive bus numbers.
  - Removed consecutive bus number requirement for `scale_load()` and
    `total_load()` functions.


Version 4.0b2 - *Mar 19, 2010*
------------------------------

#### 3/19/10
  - Released version 4.0b2.

#### 3/12/10
  - Incorporated significant updates to User's Manual (`docs/manual.pdf`).

#### 3/10/10
  - Added optional input arg to `mpver()` and other `*ver()` functions to
    trigger return of entire version struct with `Name`, `Version`,
    `Release` and `Date` (similar to MATLAB's `ver()` function).
  - Massive help text update to more closely match MathWorks conventions;
    function names in ALL CAPS, See also ..., Examples, etc.
  - Added printing of MATPOWER and MIPS version lines to verbose output.

#### 2/23/10
  - For `@opf_model`, deprecated `get_var_N()`, `get_lin_N()` and `get_nln_N()`
    methods, replaced with single `getN()` method. Added `compute_cost()`.
  - Fixed per unit bug with reserve costs and prices in `toggle_reserves()`.
  - Added `cost` field to OPF `results` struct with final values of user-defined
    costs, by named set.
  - Added `totalcost` field to `results.reserves` for OPF with reserves case,
    see `toggle_reserves()` and `runopf_w_res()`.

#### 2/2/10
  - Deprecated unused options `SPARSE_QP` and `OUT_RAW`.

#### 1/27/10
  - Renamed functions used to compute AC OPF cost, constraints and
    hessian, since they are used by more than fmincon:
     - `costfmin` --> `opf_costfcn`
     - `consfmin` --> `opf_consfcn`
     - `hessfmin` --> `opf_hessfcn`

#### 1/26/10
  - Added OPF algorithm code to output of OPF in
    `results.raw.output.alg`.

#### 1/25/10
  - Switched to using `qps_matpower()` instead of `mp_qp()`/`mp_lp()`
    for solving the DC OPF.
  - Added new top-level wrapper function for MATPOWER's QP solver,
    called `qps_matpower()`, with calling syntax similar to
    `quadprog()` from the Optimization Toolbox, to replace `mp_qp()` and
    `mp_lp()`. The main difference from the `quadprog()` API is that the
    constraints are specified as `l <= A*x <= u`, instead of
    `A*x <= b` and `Aeq*x == beq`. This new function allows for
    algorithm specific input options, return of the final objective
    function value and more detailed output reporting, such as the
    history for the trajectory returned by MIPS. The old functions,
    `mp_qp()` and `mp_lp()` are now simply wrappers around `qps_matpower()`
    and have been deprecated.
  - Added `qps_bpmpd()`, `qps_mips()` and `qps_ot()`, with interface that
    matches `qps_matpower()` to handle implementation for BPMPD_MEX,
    MIPS and Optimization Toolbox solvers, respectively.
  - Fixed a bug that could result in incorrect multipliers on
    variable bounds from the DC OPF with user-supplied linear
    constraints.

#### 1/19/10
  - Renamed the pure-MATLAB interior point solver from PDIPM to
    MIPS (MATLAB Interior Point Solver).

#### 1/18/10
  - Changed order of input args to `pdipm()`, added option for single
    input struct (like `fmincon`), more documentation, all constraints
    are now optional, returns `exitflag = -1` for 'numerically failed',
    output includes `message` field, lambda only includes relevant
    fields. Added tests for `pdipm` as standalone solver.

#### 1/12/10
  - Added saving history of trajectory of `obj`, `feascond`, `gradcond`,
    `compcond`, `costcond`, etc. for `pdipm` solver.
    See `results.raw.output.hist`.


Version 4.0b1 - *Dec 24, 2009*
------------------------------

#### 12/24/09
  - Released version 4.0b1.

#### 12/18/09
  - Make `OPF_ALG` default to 540 then 560 (no 500 MINOPF) and
    `OPF_ALG_DC` default to 200 (no 100 BPMPD_MEX).

#### 12/10/09
  - Fixed a bug, where calling `opf()` directly with individual
    input data matrices in version 2 format resulted in the matrices
    being inadvertently run through a version 1 to version 2 format
    conversion, stripping out generator capability curves, ramp
    limits and branch angle difference limits before setting up and
    running the OPF. The fix for this subtle bug involved changing
    loadcase to not assume that an input struct without a `version`
    field is in version 1 format. It now checks the size of the `gen`
    matrix to make the determination.

#### 12/8/09
  - Misc cleanup based on mlint suggestions, including:
    - Replaced `|` with `||` and `&` with `&&` where appropriate.
    - Removed unnecessary `sprintf` (and `fprintf`!?) calls from args
      to `error()`.
    - Replaced `j` (=`sqrt(-1)`) with `1j` for speed and robustness.
    - Replaced unecessary brackets `[]` with parentheses.
    - Made sure all calls to `exist()` have 2 args.
    - more

#### 12/4/09
  - Fixed bug in `savecase[]` for cases where `A` or `N` matrix is a single
    row.

#### 11/4/09
  - Removed unnecessary `return` statement at end of all M-files. If
    anything it should be an `end` statement, but even that is
    optional, so we just let functions get terminated by the
    end-of-file or another function declaration.

#### 11/3/09
  - Removed `genform.m`, `runcomp.m` and `t/t_opf.m`.
  - Renamed `compare.m` to `compare_case.m` and updated it to work with
    unsolved cases, solved PF cases and solved OPF cases.

#### 10/9/09
  - Added ability to specify interface flow limits (based on
    DC model flows).

#### 7/10/09
  - Removed `sparse_qp` and `sparse_lp` from `have_fcn()`.
  - Major speed-up in `@opf_model/linear_constraints()` for
    large problems (esp. DC OPF) and various other optimizations
    from profiling code.

#### 7/7/09
  - Fixed bug in `opf.m` introduced by automatic conversion of
    single-block piecewise linear costs to linear polynomial costs.

#### 5/27/09
  - Added `total_load.m` to provide convenient way to retreive the total
    load for the entire system, a specific zone or bus with options to
    include just fixed load, just dispatchable load, or both.

#### 5/19/09
  - The `results` struct returned by power flow or optimal power flow
    is now a strict superset of a MATPOWER case struct.
  - Extended `ext2int.m` and `int2ext.m` to handle converting entire case
    struct in a single call, storing the re-indexing info in the
    struct, using it to reorder other data structures/fields,
    execute callbacks to convert additional user data.
  - Split userfcn callbacks into multiple stages. Currently there are
    five: ext2int, formulation, int2ext, printpf, savecase.

#### 4/14/09
  - Deprecated use of `areas` data matrix. Removed it everywhere
    possible without breaking backward compatibility with version 1
    case files, which required it.
  - **INCOMPATIBLE CHANGE:** Calling `loadcase()` with 5 output arguments
    is now interpreted as ...
    - `[baseMVA, bus, gen, branch, gencost] = loadcase(casefile)`
    ... instead of ...
    - `[baseMVA, bus, gen, branch, info] = loadcase(casefile)`

#### 3/25/09
  - Added `add_userfcn.m` as to make it easy to add a new
    userfcn to a case struct, whether or not it already
    has any. Modified the fixed reserve example to use this.

#### 3/24/09
  - Added step-controlled PDIPM variant (`OPF_ALG = 565`) of
    AC OPF solver.

#### 3/19/09
  - Added `pdipm_qp()` as a new QP/LP solver based on the
    pure MATLAB PDIPM solver used for AC OPFs.
  - Added option `OPF_ALG_DC` and refactored some code to allow
    the user to select the desired solver for DC OPF.
  - Added code to `opf.m` to automatically convert single-block
    piecewise linear costs to linear polynomial costs to reduce
    the number of variables and constraints in the problem.

#### 3/17/09
  - Numerous code optimizations based on profiling code, e.g.
    changed all calls to `spdiags()` to equivalent call to `sparse()`.

#### 3/13/09
  - Added a pure MATLAB implementation of the PDIPM (primal-dual
    interior point method) solver for the AC OPF. Now the default
    solver (`OPF_ALG = 560`) if there are no optional MEX solvers
    installed.
  - Modified `fmincopf`, `copf`, `lpopf` and `dcopf` to allow branch
    `RATE_A = 0` or `RATE_A > 1e10` to mean the branch is unconstrained
    (not included at all in inequality constraints). TSPOPF solvers
    already did this. Included in tests.

#### 3/11/09
  - Allow userfcn to be an array, with elements processed in order.

#### 1/14/09
  - New version of `case39.m` with some additional versions, created
    from documented sources.

#### 7/3/08
  - Added a top level program, `runopf_w_res()`, to solve an OPF with
    fixed reserve requirements. This is a good example of how to use
    the new userfcn mechanism to add vars, costs, constraints to an
    OPF (see also `toggle_reserves.m` and `t_case30_userfcns.m`).
  - Added option to return solution as a `results` struct to `runpf()`,
    `runopf()`, `runuopf()`, `rundcpf()`, `rundcopf()`, `runduopf()`.
  - Updated `uopf.m` so input/output args match `opf.m`.
  - Added option `ENFORCE_Q_LIMS = 2` for runpf to allow one-at-a-time
    conversion of buses from PV to PQ for generator reactive power
    limit violations.
  - Fixed a (new) bug which caused the DC OPF solver to crash on
    problems with only polynomial costs.
  - Added `userdata` to `@opf_model` object.

#### 6/10/08
  - Added new way to specify user vars, constraints, costs via
    userfcn for OPF.
  - Added option to return OPF `results` in a struct.
  - Added defaults for user cost params in `fparm` and `H`, making them
    optional even when `N` and `Cw` are given.

#### 5/22/08
  - Major refactorization of OPF implementation with shared code
    for a generalized formulation that includes the DC OPF as
    well as the legacy solvers based on `constr` and `LPconstr`.
  - Deprecated `OPF_ALG` values 100, 120, 140, 160, 200, 220, 240,
    and 260 in favor of the new generalized formulation
    equivalents 300, 320, 340 and 360.
  - Removed options `OPF_ALG_POLY`, `OPF_ALG_PWL` and
    `OPF_POLY2PWL_PTS`.

#### 5/2/08
  - Move OPF input argument processing into `opf_args.m`, now
    shared by `opf.m`, `dcopf.m`, `fmincopf.m`, `mopf.m` and `tspopf.m`.
  - Rewrote the DC OPF to include generalized user constraints,
    costs and extra vars (like AC formulation). Note, that if
    `A` or `N` have enough columns for the AC formulation, `opf.m`
    assumes they are for the AC OPF and strips out the extra
    columns before passing to `dcopf.m`.
  - Added the ability to read and save generalized OPF user
    constraints, costs and var limits in case struct.
  - Modified `savecase.m` to include saving of `MU_ANGMIN`, `MU_ANGMAX`
    columns of `branch` matrix.

#### 3/13/08
  - Added a function `makeLODF.m` to compute line outage distribution
    factors.
  - Added a function `scale_load.m` to scale load by zones.

#### 3/7/08
  - Updated `fmincopf` and `mpoption` to work with version 4 of
    Optimization Toolbox. Added option `FMC_ALG` for select between
    fmincon's active set method and variations of the new
    interior-point algorithm.
  - Added functions to compute second derivatives of constraints
    and cost (explicit Hessian evaluation) for use with
    interior-point solvers, etc.
  - **INCOMPATIBLE CHANGE:** `dAbr_dV()` now gives partial derivatives
    of the *squared* magnitudes of flows w.r.t. V, as opposed
    to the magnitudes.
  - Modified the implementation of all flow constraints for `fmincon`
    (and `constr`) to use squared flow limits instead of absolute
    value to correctly avoid div-by-zero errors in computing
    gradients, and to prepare for implementing Hessian code.
    Shadow prices still correspond to absolute value limits.
  - Fixed bug in `fmincon` (and `constr` and LP) based OPF which
    allowed an active power flow limit to be violated when using
    `OPF_FLOW_LIM = 1` (missing absolute value).

#### 3/3/08
  - **INCOMPATIBLE CHANGE:** Changed input argument order for `uopf`
    and added general linear constraints and generalized costs.

#### 1/10/08
  - Significant speed improvements in `makeYbus.m` and `makeBdc.m`.


Version 3.2 - *Sep 21, 2007*
----------------------------

#### 9/21/07
  - Released version 3.2.

#### 9/17/07
  - Added option to `cdf2matp.m` to specify output case file version.

#### 9/7/07
  - Fixed bug in `pfsoln.m` which caused incorrect value for `Qg` when
    `Qmin == Qmax` for all generators at a bus in power flow solution.
  - Added 5 larger scale (> 2000 bus) cases for Polish system.
    *Thanks to Roman Korab <roman.korab@polsl.pl>.*
  - Modified default OPF algorithm selection to use PDIPMOPF
    if available and MINOPF is not. Order of precedence is now
    500, 540, 520, 100/200.

#### 7/6/07
  - Added ability in `opf.m` and `fmincopf.m` to specify initial value
    and bounds on user variables via new input arguments `z0`, `zl`, `zh`.

#### 6/22/07
  - **INCOMPATIBLE CHANGE:** Name of option 24 in mpoption change from
    `OPF_P_LINE_LIM` to `OPF_FLOW_LIM`.
  - Added option to use current magnitude instead of apparent power
    for line flow limits. Set `OPF_FLOW_LIM` to 2.

#### 6/21/07
  - **INCOMPATIBLE CHANGE:** Changed the sign convention used for
    phase shifters to be consistent with PTI, PowerWorld, PSAT, etc.
    E.g. A phase shift of 10 deg now means the voltage at the "to"
    end is delayed by 10 degrees.

#### 6/15/07
  - Added `t_auction_pdipm.m` and renamed `t_auction.m` to
    `t_auction_minopf.m`.

#### 6/8/07
  - Updated `have_fcn.m` to check for appropriate minimum versions of
    MATLAB, for TSPOPF.

#### 6/7/07
  - Modified `printpf.m` to correctly detect binding line limits when
    a limit of 0 is taken to mean unconstrained.
  - Fixed bugs in handling of multipliers for general PQ capability
    curves in `fmincopf.m` (also in `mopf.m` and `tspopf.m`).
  - Refactored `t_opf.m` into separate files for each solver.
  - Modified `opf.m`, `mpoption.m`, `mpver.m`, `have_fcn.m` to include
    support for TSPOPF, a new optional package of OPF solvers.

#### 9/29/06
  - Added check to `runpf.m` for case where all gens hit Q limits when
    `ENFORCE_Q_LIMS` is enabled.


Version 3.1b2 - *Sep 15, 2006*
------------------------------

#### 9/15/06
  - Released version 3.1b2.

#### 9/12/06
  - Added `makePDFT.m` which builds the DC PTDF matrix for a specified
    slack distribution.

#### 8/16/06
  - Added optional outputs `xr`, `pimul` to `fmincopf` and `opf.m` to make them
    fully interchangeable with `mopf.m`.

#### 8/15/06
  - Added branch angle difference constraints to general OPF formulation
    in `fmincopf.m` (and `mopf.m`). These limits are specified by non-zero
    values in the `ANGMIN` and/or `ANGMAX` columns of the `branch` matrix.
    If limits are provided in the data, they are enforced by default.
    This can be overridden by setting the `OPF_IGNORE_ANG_LIM` option
    to 1 using `mpoption`.
  - Fixed (invisible) bug with multipliers of lower bounded linear
    constraints in `fmincopf.m`.


Version 3.1b1 - *Aug 1, 2006*
-----------------------------

#### 8/1/06
  - Released version 3.1b1.

#### 4/28/06
  - Fixed `mpver.m` so it will properly handle case where the Optimization
    Toolbox is not installed.

#### 3/15/06
  - **INCOMPATIBLE CHANGE:** Updated `opf.m`, `fmincopf.m`, `costfmin.m`, `consfmin.m` to
    be able to be compatible with latest MINOPF. User supplied A matrix for
    general linear constraints no longer includes columns for `y` variables
    (helper vars for piecewise linear gen costs), and now requires columns
    for all `x` (OPF) variables. Added generalized cost model and generator PQ
    capability curves.
  - Modified `savecase.m` to always save MAT files with -V6 under newer MATLAB
    versions.
  - Added a number of tests to `t_opf.m` for MINOPF and fmincopf for generalized
    costs and additional linear constraints. Added test for fmincopf for
    generator PQ capability curves.

#### 3/10/06
  - Added baseKV data to `case118.m` from PSAP file 
    <http://www.ee.washington.edu/research/pstca/pf118/ieee118psp.txt>.

#### 3/8/06
 - Renamed col 5 of `gencost` from `N` to `NCOST` everywhere.

#### 10/14/05
  - Updated version 2 case file format to modify generator PQ capability
    curve specifications.
  - Added `hasPQcap.m` and test for gen PQ capability curve in OPF.

#### 8/22/05
  - Added `OPF_IGNORE_ANG_LIM` option to `mpoption.m`.

#### 8/5/05
  - Modified identification of binding constraints in `printpf.m`. A
    constraint is now considered to be binding if the tolerance is less
    than or equal to `OPF_VIOLATION` tolerance -OR- if the corresponding
    Kuhn-Tucker multiplier is non-zero. This allows binding generator
    capability curves to be reported via multipliers on Pg and Qg limits.

#### 7/8/05
  - Updated `loadcase.m`, `savecase.m`, `idx_bus.m`, `idx_gen.m`, `caseformat.m`
    and tests for version 2 case file format, which includes piece-wise
    linear generator capability curves, generator ramp rates and branch
    angle difference limits.
  

Version 3.0.0 - *Feb 14, 2005*
------------------------------

#### 2/14/05
  - Released version 3.0.0.

#### 2/3/05
  - In `mp_lp.m` and `mp_qp.m`, on Windows it now makes sure BPMPD_MEX is not
    called in verbose mode which causes a MATLAB crash.


Version 3.0b4 - *Jan 28, 2005*
------------------------------

#### 1/28/05
  - Released version 3.0b4.

#### 1/27/05
  - Added `case6ww.m` and `case4gs.m`.
  - Minor modifications to `printpf.m` to handle larger bus numbers.

#### 1/26/05
  - Minor changes to `uopf.m` to make sure it plays nicely with dispatchable
    loads.

#### 1/25/05
  - Major updates to user manual.

#### 1/24/05
  - Switched to using the new `isload()` to check for dispatchable load.
  - For dispatchable loads, switched from using `PG` and `QG` to `PMIN` and either
    `QMIN` (for inductive loads) or `QMAX` (for capacitive loads) to define the
    constant power factor constraint. This prevents the power factor
    information from being lost when it is dispatched to zero. If the initial
    values of `PG` and `QG` are not consistent with the ratio defined by `PMIN`
    and the appropriate Q limit it gives an error. This is to prevent a user
    from unknowingly using a case file which would have defined a different
    power factor constraint under previous versions of MATPOWER.
    If both `QMIN` and `QMAX` are zero, it no longer includes the redundant
    unity power factor constraint.

#### 1/20/05
  - Updated `printpf.m` to display dispatchable loads and generators
    separately. Reorganized the area summary section and corrected the net
    exports value (subtracted half of tie-line loss) to make the numbers
    add up correctly.

#### 1/18/05
  - Added to `runpf.m` the ability to enforce generator reactive power limits
    by allowing the voltage to deviate from the set-point. This option is
    controlled by the new `ENFORCE_Q_LIMS` option, which is off by default.
    *Thanks to Mu Lin of Lincoln University, New Zealand whose contributions
    inspired this feature.*
  - Modified `pfsoln.m` to divide reactive power dispatch between multiple
    generators at a bus in proportion to each gen's reactive power range,
    as opposed to equally. This means that all generators at a bus will
    reach their upper (or lower) limits simultaneously.
  - Added generator status column to generator section of `printpf.m` output.
    Fixed bugs where non-zero output of decommitted generators was displayed
    and included in generation totals in generator and bus sections.

#### 1/14/05
  - Moved some setting of `MNS_*` default options from `opf.m` to `mopf.m`.
  - Eliminated unused output args in `dcopf.m`.
  - Modified `printpf.m` to zero out reactive generator output for DC cases
    and to use `OPF_VIOLATION` tolerance to detect binding constraints, as
    opposed to non-zero Kuhn-Tucker multipliers.

#### 1/12/05
  - Modified bpmpd portion of `mp_qp.m` and `mp_lp.m` to use default value for
    `TFEAS2` and eliminate variable limits which appear to be artificial
    large values used to indicate free variables.

#### 1/4/05
  - Fixed potential bug in dimensions of `Yf` and `Yt` created in `makeYbus.m`.

#### 12/17/04
  - Added feasibility check to `mp_lp.m` and `mp_qp.m` to work around a
    recently discovered bug in BPMPD_MEX 2.21 where it sometimes returns an
    incorrect (infeasible) solution for a DC OPF problem. This bug has yet
    to be encountered in any other context.

#### 12/13/04
  - Added `mpver.m` to print version information.

#### 9/23/04
  - Fixed bugs in `cdf2matp.m` which prevented it from working at all
    when not specifying both input parameters and caused it to
    sometimes not add the warnings at the end of the file.
  - Fixed typo in name of lower bound input argument in `opf.m`. Only
    affected those calling OPF directly with extra linear constraints.


Version 3.0b3 - *Sep 20, 2004*
------------------------------

#### 9/20/04
  - Released version 3.0b3.
  - Generated clean versions of all included case files using latest
    `cdf2matp` and `savecase`. Added documentation for source of data
    for case files.
  - More enhancements to `cdf2matp.m`. Adds comments at beginning, appends
    conversion warnings as comments at end of file. Uses `savecase.m` to
    save the data.
  - Updated `savecase.m` to use `%g` instead of `%f` many places, correctly
    handle multi-line comments, include headers for extra columns for
    solved cases. Optionally returns filename with extension.

#### 9/17/04
  - Fixed bug in `grad_std.m`, introduced in 3.0b2, which prevented `constr`
    and LP-based OPF solvers from working for polynomial cost functions.

#### 9/15/04
  - In `cdf2matp.m`, added input args, updated docs, switched to named
    indexing of data matrices, new method for creating gen costs.
  - Documentation fixes and additions from Pan Wei.


Version 3.0b2 - *Sep 7, 2004*
-----------------------------

#### 9/7/04
  - Released version 3.0b2.
  - Added `OPF_P_LINE_LIM` option to mpoptions to use active power
    instead of apparent power for line limits. *Thanks to Pan Wei
    for the suggestion and some code.*

#### 9/1/04
  - Fixed bug in `savecase.m` introduced when making `areas` and `gencost`
    optional.
  - Updated `opf_slvr.m` with options for MINOS and fmincon.
  - Removed option 15 `OPF_NEQ` from docs (not a user option). Removed option
    52 `VAR_LOAD_PF` (unused, always behaves as if this option were 1).
    Changed semantics and default value for option 51 `SPARSE_QP`. By default
    (value = 1) it will use sparse matrices if a sparse QP/LP solver is
    available and full matrices otherwise. Setting the value to 0
    will force it to use full matrices even with a sparse-capable solver.
  - Cleaned up checking for optional functionality, and fixed a bug
    that would miss MEX files if there was an identically named directory
    by adding `have_fcn.m`.


Version 3.0b1 - *Aug 25, 2004*
------------------------------

#### 8/25/04
  - Released version 3.0b1.

#### 8/24/04
  - Made `mpoption()` throw an error if passed an invalid option name.

#### 8/23/04
  - Added an `fmincon` based OPF solver for the generalized formulation
    previously used by `mopf` (Carlos).
  - Restructured `opf.m` so all OPF solvers have a similar API based
    on the one from `mopf.m` (Carlos).
  - Added some quick tests for `runpf` and `runopf` for each algorithm.

#### 8/13/04
  - Renamed `area` variable to `areas` to avoid masking the built-in
    function of the same name.
  - Made OPF data matrices `areas` and `gencost` optional for running
    simple power flow.

#### 7/15/04
  - The `loadcase` function (and therefore all of the `run*` functions
    now optionally accept a struct with the data matrices as fields
    in place of the case file name.
  - Added `t` subdirectory with various tests and testing tools.

#### 7/8/04
  - Updated `mp_lp.m` and `mp_qp.m` to try `linprog()` and `quadprog()`
    after trying `bp`, since `lp()` and `qp()` are no longer included
    in the Optimization Toolbox as of version 3.

#### 7/7/04
  - Removed `case.m`, added `caseformat.m`, made `case9.m` the default
    case and fixed function names in other case files to avoid
    use of reserved word `case`.
  - Fixed bugs in `runcomp.m`.

#### 6/23/04
  - Fixed bug in `newtonpf.m` which caused algorithm to diverge when
    the Newton step resulted in a negative voltage magnitude.

#### 4/17/03
  - Changed `uopf.m` to use a dynamic programming approach. More
    computationally expensive, but should find significantly better
    results when there are many gens to shut down.
  - Added `mp_lp.m` and `mp_qp.m`, equivalents to `lp.m` and `qp.m`,
    respectively that call BPMPD if available. Modified `LPrelax.m`,
    `LPsetup.m` and `dcopf.m` to call these new functions.

#### 4/14/03
  - Fixed a bug in `pfsoln.m` which for cases with a single generator.

#### 10/23/02
  - Fixed bus numbering bug in System Summary section of `printpf.m`.

#### 6/20/00
  - Fixed a bug in `printpf.m` in the generator section, where
    the generator was assumed to be off if it's real power
    output was zero, even if the reactive output was non-zero.
  - Modified `printpf.m` to print out lambdas in generation section
    for generators that are shut down.

#### 6/8/00
  - Modified `cdf2matp.m` so that `Pd` also includes any generation at
    buses specified as PQ buses. Also modified identification of
    generator buses to include only PV or reference buses. *Thanks
    to Venkat.*
  - Modified `cdf2matp.m` so that it always treats the input values
    for `Gs` and `Bs` as per unit values and converts them to actual
    values expected by MATPOWER. *Thanks to D. Devaraj.*

#### 11/9/99 - version 2.5b3

#### 9/22/99
  - Modified `grad_*.m` to return sparse matrices, unless using
    `constr.m` or an LP/QP solver that doesn't handle sparse
    matrices. Cleaned up sparse<->full conversions in `LPconstr.m`,
    `LPrelax.m`, and `LPsetup.m`.

#### 9/21/99
  - Undid a "bug fix" from 3/6/98 in `makeYbus.m` which zeros out
    charging capacitance for transformers. Apparently some
    transformer models actually have a non-zero charging parameter
    when using the model used by MATPOWER (ideal transformer in
    series with a PI model).
  - Added `loadcase.m` which loads a MATPOWER case from an M-file
    or from a MAT-file. Changed all of the `run*.m` files to use this
    as the default way to load case files.
  - Renamed `print2mp.m` to `savecase.m` and added the ability to
    save a case as a MAT-file as well as an M-file.

#### 9/15/99
  - Fixed `opf.m` so that it correctly uses the termination
    tolerances in the MATPOWER options vector for `constr.m`.
  - In previous versions, Pmin/Pmax constraints are relaxed by
    `10 * OPF_VIOLATION` in `opf.m` to avoid falsely reporting a
    binding Pmin/Pmax constraint in a case where a piece-wise linear
    cost function has a corner point exactly at the limit. This
    code was moved out of `opf.m` (and the standard MATPOWER
    distribution) to `smartmkt.m` and the value was changed to
    `100 * OPF_VIOLATION`.
  - Modified `opf.m` so the MINOS-based solver uses `OPF_VIOLATION`
    to set the value of `MNS_FEASTOL` and `MNS_ROWTOL` if they are
    set to zero.

#### 9/9/99
  - Included MINOS-based OPF with all of its options as
    algorithm 500. (involved including `area` in calls to `opf.m`
    and `uopf.m`)
  - Removed some unused lines from `fun_ccv.m` and `grad_ccv.m`.

#### 8/5/99
  - Fixed a bug in the `pfsoln.m` in the distribution of Q among
    generators at the same bus. Initially attempted to distribute
    Q to generators proportional to each generators' Q "capacity".
    To do this correctly requires special cases for generators
    with `QMIN` equal to `QMAX`. For the sake of simplicity, we now
    distribute Q equally among all generators at the bus.
    Note: As before, the simple power flow does NO feasibility
    checking.

#### 7/19/99
  - Modified `runuopf.m` and `uopf.m` to handle DC OPF. Added the
    function `runduopf.m` which calls `runuopf.m` with the `PF_DC` flag
    set to 1.
  - Fixed size of 2nd order (all zero) coefficient of objective
    for piecewise linear cost case in `dcopf.m`.

#### 7/16/99
  - Added the flag `QP_SPARSE` to `mpoption.m` to indicate whether the
    QP solver being used can accept sparse matrices. Also modified
    `dcopf.m` to use this flag.
  - Fixed handling of `VERBOSE` option in `dcopf.m`
  - Added the flag `PF_DC` to `mpoption.m` to indicate whether the
    power flow formulation to be used for power flow and optimal
    power flow is a DC approximation or full AC representation.
    Merged `rundcpf.m` with `runpf.m` and `rundcopf.m` with `runopf.m`
    so that the appropriate solver will be used based on the
    value of the `PF_DC` flag in the options. The functions `rundcpf.m`
    and `rundcopf.m` were modified to simply call `runpf.m` and
    `runopf.m`, respectively, with the `PF_DC` flag set to 1.

#### 7/15/99
  - Changed the sign of the phase shifters in `printpf.m` to be
    consistent with the bug fix to `makeYbus.m` made on 3/6/98.

#### 7/14/99
  - Included four new m-files (`makeBdc.m`, `dcopf.m`, `rundcpf.m`,
    and `rundcopf.m`) which implement a DC power flow and DC
    optimal power flow algorithms.

#### 7/13/99
  - Cleaned up variable names in `makeYbus.m` to avoid confusion.

#### 6/10/99
  - Changed UOFP to UOPF in print statements `uopf.m`.

#### 6/3/99
  - Modified `print2mp.m` overwrite instead of append to an
    existing file.
  - Fixed bug in `cdf2matp.m` to make it always correctly write
    a text file output.

#### 6/2/99 - version 2.5b2
  - Modified `print2mp.m` to include line flows and Lagrange
    and Kuhn-Tucker multipliers in output if available.

#### 4/29/99
  - Included a Gauss-Seidel power flow solver `gausspf.m`, and
    made corresponding changes to `runpf.m` and `mpoption.m`.
    *Code contributed by Alberto Borghetti.*

#### 4/28/99
  - Modified `newtonpf.m` to handle cases with no PQ buses or no
    PV buses under newer versions of MATLAB.

#### 2/25/99
  - Fixed a bug in `uopf.m` which occurs when two (or more)
    generators have positive decommitment indices but shutting
    them down one at a time always results in increased system
    cost. In this scenario, it would go into an infinite loop
    of attempting to shut them down one by one.

#### 2/24/99
  - Modified `uopf.m` to be able to handle the case where the
    sum of the Pmin's is greater than the load. It shuts down
    generators in order of decreasing average cost at Pmin
    (breaking ties randomly) until this infeasibility is gone.

#### 2/16/99
  - Fixed bug in `pfsoln.m` which caused crashes in MATLAB 5
    for systems with no capacitors.
  - Added `print2mp.m`, which can print out a MATPOWER case file
    from the data matrices.
  - Added to run*`pf.m` ability to save solved case.

#### 2/10/99
  - Modified `ext2int.m` to allow for area matrix to be empty.

#### 12/3/98
  - Changed `pfsoln.m` so that there is only one slack generator.
    Instead of dividing the P between multiple gens at the
    slack bus in proportion to capacity (this caused problems
    for the `LPconstr` versions of the OPF), it now treats the
    first generator at the slack bus as the only slack generator,
    leaving the dispatch of the other gens at the bus unchanged.
  - Added generator number to generation constraint printout and
    branch number to branch data and branch flow limit printouts.

#### 12/2/98
  - Changed `printpf.m` to print elapsed time and objective fcn
    value even when `OUT_SYS_SUM` is turned off.
  - Added code to `LPconstr.m` to explicitly zero out lambdas for
    non-binding constraints.

#### 12/1/98
  - Made modifications to the following to allow for multiple
    generators at each bus. For simple power flow, the Q dispatch
    is divided between multiple gens at a bus in proportion to
    each gen's Q capacity. Likewise with P for multiple gens at
    the slack bus.
    - `bustypes.m`
    - `fun_ccv.m`
    - `fun_std.m`
    - `grad_ccv.m`
    - `grad_std.m`
    - `LPeqslvr.m`
    - `makeSbus.m`
    - `opf.m`
    - `opfsoln.m`
    - `pfsoln.m`
    - `printpf.m`
    - `runpf.m`

#### 10/29/98
  - Fixed bug in `uopf.m` which caused it to crash when attempting
    to restart a generator after more than 2 had been shut down.

#### 10/19/98
  - Generalized definition of `GEN_STATUS` column of `gen` matrix
    to allow for distinctions in the status of out-of-service
    generators. The default values of 0 => out-of-service and
    1 => in-service still work, but the logic has been changed
    so that `GEN_STATUS > 0` is now in-service and
    `GEN_STATUS <= 0` is now out-of-service, as opposed to
    `GEN_STATUS ~= 0` and `GEN_STATUS == 0`, respectively, which
    was used previously. This allows for a `GEN_STATUS` of -1,
    for example, to indicate a generator which is off-line
    but could be brought on in case of an emergency.

#### 9/2/98
  - Fixed bug in `printpf.m` which caused area exports to be
    off slightly.

#### 9/1/98
  - Fixed bug in `printpf.m`. Total intertie flow was double the
    correct value.

#### 8/31/98
  - Fixed bug which included line flow limits for out-of-service
    lines in OPF.
  - Modified `pfsoln.m`, `opfsoln.m`, `printpf.m` to zero out flow on
    lines which are out-of-service. *Found by Ramazan Caglar.*

#### 7/28/98
  - Changed VAR and MVAR to VAr and MVAr everywhere in output.

#### 3/13/98
  - Decreased the default value of `LPC_TOL_X` option to increase
    solution quality.
  - Modified fix of 2/10/98 to use a value based on the value of
    the `OPF_VIOLATION` option.

#### 3/6/98
  - Fixed 2 bugs in `makeYbus.m`. Phase shifters now shift the phase the
    right direction, the line charging susceptance parameter is now
    correctly ignored for transformer and phase shifters.

#### 3/3/98
  - Fixed a bug `fun_std.m` which caused it to always compute 2nd order
    derivatives. Now it only computes them when requested.

#### 2/10/98
  - In previous versions, Pmin/Pmax constraints are relaxed by 1.0e-6
    in `opf.m` to avoid falsely reporting a binding Pmin/Pmax constraint
    in a case where a piece-wise linear cost function has a corner
    point exactly at the limit. Changed the amount of relaxation to
    1.0e-4 since the problem still occurred at times.

#### 1/29/98
  - Changed the value of `LPC_MAX_IT` from 1000 to 400 to allow for
    earlier detection of infeasible OPF.

Version 2.0 - *Dec 24, 1997*
----------------------------

#### 12/24/97
  - Released version 2.0.

#### 12/19/97
  - Fixed ambiguity in case file data and comments regarding lines
    vs. transformers. Now a tap ratio of zero means that it's a line
    and a non-zero tap ratio means that it's a transformer.
  - Optimized formation of Ybus (and hence B matrices).

#### 12/18/97
  - Implemented fast decoupled load flow.

#### 12/17/97
  - Optimized formation of Jacobian matrix in `newtonpf.m` (significant
    improvement for large systems under MATLAB 5).
  
#### 12/16/97
  - Fixed another bug in calculation of losses. Previous versions
    did not take into account off-nominal taps for transformers.
  - Fixed a bug in calculation of losses. Previous versions
    included line charging injection in reactive line losses.
  - Added ability to optionally return solution data from
    `run*.m` functions.
  - Added ability to optionally print results to a file.
  - Added system and area summaries to `printpf` and modified to
    handle the new printing options.

#### 12/12/97
  - Consolidated printing into `printpf.m`, eliminated `printopf.m`.
  - Removed QCCV method (standard formulation solves same problem,
    but more efficiently).
  - Removed OPF algorithms which use fixed generator voltages
    (this can still be done by changing voltage limits in the
    case file), renumbered OPF algorithms, removed `CCV.m` and
    `varVg.m`.

#### 12/11/97
  - Added 2 more levels of control of verbose output.
  - Put all MATPOWER options into an options vector defined in
    `mpoption.m`.

#### 12/10/97
  - Incorporated new LP-based OPF routines and updated alg codes.
  - Fixed a bug in the documentation in the case files regarding
    the 4th column of `gencost`. For piece-wise linear cost functions
    this value is the number of data points, not the number of
    parameters (x and y for each point).
  - Removed some m-files that are not used (`usesOT.m`, `usesLP.m`).
  - Renamed some m-files (`OTfungra.m` to `fg_names.m`, `OTSfun.m` to
    `fun_std.m`, `OTgra.m` to `grad_std.m`, `OTCCVfun.m` to `fun_ccv.m`,
    `OTCCVgra.m` to `grad_ccv.m`).

#### 12/8/97
  - Rewrote `uopf.m` to use a smarter decommitment strategy (see the
    docs for the details of the new method). Removed `ref`, `pv`, `pq`
    from the list of parameters passed in, since they were not used.

#### 11/19/97
  - Fixed a bug in previous versions of `uopf.m` which returned
    incorrect values for Pmin.

#### 10/28/97
  - Increased maximum number of iterations for constr-based OPF.

#### 10/24/97
  - Fixed a bug in previous versions which may result in incorrectly
    reporting Pmin or Pmax limits to be binding, possibly with large
    multipliers, if the piece-wise linear cost function has a corner
    point exactly at Pmin or Pmax.

#### 10/22/97
  - Added to `OTSgra.m` (renamed to `grad_std.m` in 2.0) the ability
    to return the second derivatives of the objective function.

#### 9/24/97
  - Fixed a bug in previous versions of `runuopf.m` which prevented it
    from printing out the raw data needed for our Perl DB interface.

#### 9/23/97
  - Fixed a bug in 1.1b1 in `OTCCVgra.m` (renamed to `grad_ccv.m` in 2.0)
    which caused printing of warning message "Concatenation involves
    an incommensurate empty array" under MATLAB 5.

#### 9/22/97
  - Fixed a bug in 1.1b1 which prevented `runuopf.m` from running at all.
    Wrong number of parameters to call `opf.m`.

#### 9/20/97
  - Released version 1.1b1.

#### 9/19/97
  - Modified the formulation of the OT-based OPF. The objective
    function may now include costs for reactive power as well as
    active power. In previous versions the reactive power variables
    and reactive power balance equations for generator buses were
    not included explicitly in the optimization as variables and
    equality constraints. Generator reactive powers were computed
    directly. Now they are included explicitly in the optimization.
    Costs for Qg are specified in extra rows in `gencost`.


Version 1.0.1 - *Sep 20, 1997*
------------------------------

#### 9/20/97
  - Released version 1.0.1.

#### 9/19/97
  - Fixed a bug in 1.0 `OTSgra.m` and `OTCCVgra.m` (renamed to
    `grad_std.m` and `grad_ccv.m`, respectively, in 2.0) which used
    incorrect coefficients to compute cost if specified as
    polynomials of different degrees.

#### 9/18/97
  - Fixed a bug in 1.0 in `OTopf.m` which caused the last equality
    constraint (Q mismatch for last PQ bus) to be treated as an
    inequality constraint. It appears that this constraint was
    normally binding (unless Qd was negative) in which case the
    solution was still correct.
  - Fixed a bug in 1.0 in `runpf.m`, initial voltage for generators
    which were shut down were taken from `gen(:, VG)` rather
    than `bus(:, VM)`.
  - Fixed a bug in 1.0 in `varVg.m` which caused Kuhn-Tucker
    multipliers to print out in the wrong place for LP-based OPF.


Version 1.0 – *Sep 17, 1997*
----------------------------

#### 9/17/97
  - Released version 1.0 (first widely publicized release).
  - added placeholders for LP-solvers that we can't re-distribute
  - updated documentation

#### 9/12/97
  - added ability to do pretty & ugly printing at the same time
    also documented that ugly printing is for talking to our
    our Perl database interface code
  - included Deqiang (David) Gan's LP IEEE -> matpower data
    conversion code
  - included Deqiang (David) Gan's LP based OPF code
  - fixed `LAM_Q` bug, now computes correctly for generator buses
  - fixed some bugs in `totcost.m`

#### 9/9/97
  - removed `PRICE` from `idx_gen`

#### 9/4/97
  - added code to convert from (possibly non-consecutive) external
    bus numbering to consecutive internal bus numbering before
    solving, and back to external before printing results
  - replaced `test*pf` with `run*pf` which are now functions
    taking the casefile name as a parameter (among other params)
  - made changes necessary to handle new format of case file
    (generator costs moved to `gencost` variable)


First Public Release – *Jun 25, 1997*
-------------------------------------

#### 6/25/97
  - made first public release (not widely publicized)
  - documentation updates
  - changed names of m-files to fit DOS 8.3 limitation
    - `buildsbus.m`     =>  `makeSbus.m`
    - `buildybus.m`     =>  `makeYbus.m`
    - `idx_branch.m`    =>  `idx_brch.m`
    - `dSbranch_dV.m`   =>  `dSbr_dV.m`
    - `dAbranch_dV.m`   =>  `dAbr_dV.m`
    - `ucopfsoln.m`     =>  `uopfsoln.m`
    - `testucopf.m`     =>  `testuopf.m`
    - `ucopf.m`         =>  `uopf.m`  (for naming consistency)
  - changed copyright notice

#### 6/18/97
  - modified `ucopf.m` to allow a generator to be turned back on if
    shutting it off results in an infeasible (or at least
    non-convergent) OPF, also changed the order of shutting down
    generators which are dispatched at zero, now chooses one with
    largest `mu_Pmin`

#### 6/12/97
  - fixed bug in `printpf.m` so it doesn't print `PG` & `QG` for gens that
    have been shut down
  - fixed bug in `pfsoln.m` to correctly compute the reference bus power
    injection when generators have been shut down

#### 6/10/97
  - fixed Vg initialization bug in `testpf.m` (not just `testopf`, etc)

#### 6/9/97
  - fixed bug in PLCCV versions which set the initial values of the
    cost variables wrong (used p.u. Pg instead of actual)
  - made `opfsoln.m` copy generator voltages back to `gen(:, VG)`
  - fixed bug in code which initializes generator voltages, it was
    always setting the angle to zero, now it uses the value from the
    case file

#### 6/3/97
  - included OPF variations which use cost variables constrained
    by a piece-wise linear cost function (PLCCV = piece-wise linearly
    constrained cost variables)

#### 6/2/97
  - included OPF variations which use cost variables constrained
    by a quadratic cost function (QCCV = quadratically constrained
    cost variables)
  - included OPF variation which allows generator voltage
    magnitudes to vary
  - fixed line in `test*pf.m` scripts which initializes `V0` (I'd missed
    the `sqrt(-1)` before

#### 4/16/97
  - changed line 59 of `ucopf.m` from `return` to `break` to ensure
    return values are correct

#### 4/14/97
  - added some print statements to `ucopf.m`

#### 4/12/97
  - reduced max iterations to 100 for `constr` in `opf.m`

#### 4/8/97
  - modified `opf.m`, `ucopf.m`, `testopf.m`, `testucopf.m` to include
    `success`, a variable which indicates whether OPF was solved
    successfully or not

#### 4/7/97
  - fixed bug in `ucopf.m`, assumed all generators are initially
    available

--

[1]: https://github.com/MATPOWER/mptest
[2]: https://github.com/MATPOWER/mips
[3]: https://github.com/MATPOWER/most
