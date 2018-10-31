What's New in MATPOWER 7.0
--------------------------

#### Released Oct 31, 2018

Below is a summary of the changes since version 6.0 of MATPOWER. See the
[`CHANGES.md`][1] file for all the gory details. For release notes for
previous versions, see Appendix H of the [MATPOWER User's Manual][2].

#### New Features:
- New MATPOWER installer script `install_matpower()` automatically
  updates MATLAB or Octave paths or optionally provides the commands
  required to so.
- Support for additional user-defined general nonlinear constraints and
  costs in AC OPF.
- Support for exporting MATPOWER case to PSS/E RAW data format.
- Three new variants of the standard AC OPF formulation, for a total of
  four, including both nodal power and current balance constraints and
  both polar and cartesian representations of voltage. See the new
  `opf.current_balance` and `opf.v_cartesian` options.
  *Thanks to Baljinnyam Sereeter.*
- Three new power flow algorithms for radial distribution systems
  selected via the three new options for `pf.alg`, namely `'PQSUM'`,
  `'ISUM'`, `'YSUM'`. Also includes new MATPOWER options
  `pf.radial.max_it` and `pf.radial.vcorr`. See Section 4.3 on
  "Distribution Power Flow" in the [MATPOWER User's Manual][2] for
  details. *Thanks to Mirko Todorovski.*
- Major update to OPF soft limit functionality, supporting soft limits
  on all AC and DC OPF inequality constraints, including branch flow
  constraints, bus voltage bounds, generator active and reactive
  bounds, branch flow and branch angle difference limits.
  *Thanks to Eran Schweitzer.*
- New options:
  - `pf.nr.lin_solver` controls the linear solver used to compute the
    Newton update step in the Newton-Raphson power flow.
  - `pf.radial.max_it` and `pf.radial.vcorr` are options for the new
    radial power flow algorithms.  *Thanks to Mirko Todorovski.*
  - `cpf.enforce_flow_lims` and `cpf.enforce_v_lims` control
    enforcement of branch flow and bus voltage magnitude limits in the
    continuation power flow and `cpf.flow_lims_tol` and `cpf.v_lims_tol`
    control the respective detection tolerances.
    *Thanks to Ahmad Sadiq Abubakar and Shrirang Abhyankar.*
  - `opf.current_balance` and `opf.v_cartesian` control formulation
    used for AC OPF. *Thanks to Baljinnyam Sereeter.*
  - `opf.softlims.default` determines whether or not to include soft
    limits on constraints whose parameters are not specified explicitly
    in the `mpc.softlims` struct. For use with enhanced
    `toggle_softlims()` functionality. *Thanks to Eran Schweitzer.*
  - `opf.start` replaces deprecated `opf.init_from_mpc` and adds a new
    possibility to automatically run a power flow to initialize the
    starting state for the OPF.
- New functions:
    - `calc_branch_angle` calcultes voltage angle differences for
      branches.
    - `dImis_dV()` evaluates the partial derivatives of nodal current
      balance with respect to bus voltages.
    - `d2Imis_dV2()` evaluates the 2nd derivatives of nodal current
      balance with respect to bus voltages.
    - `d2Imis_dVdSg()` evaluates the 2nd derivatives of nodal current
      balance with respect to bus voltages and generator injections.
    - `d2Abr_dV2()` evaluates the 2nd derivatives of squared branch
      flows with respect to bus voltages.
    - `gentypes()` and `genfuels` provide list of standard generator
      unit types and fuel types, respectively.
    - `install_matpower()` installs MATPOWER by automatically modifying
      or, alternatively, showing needed modifications to MATLAB or
      Octave path.
    - `loadshed()` computes MW curtailments of dispatchable loads.
    - `opf_branch_ang_fcn()` evaluates AC branch flow limit constraints
      and gradients.
    - `opf_branch_ang_hess()` evaluates Hessian of AC branch flow limit
      constraints.
    - `opf_current_balance_fcn()` evaluates AC current balance
      constraints and gradients.
    - `opf_current_balance_hess()` evaluates Hessian of AC current
      balance constraints.
    - `opf_veq_fcn()` evaluates voltage magnitude equality
      constraints and gradients.
    - `opf_veq_hess()` evaluates Hessian of voltage magnitude equality
      constraints.
    - `opf_vlim_fcn()` evaluates voltage magnitude limit constraints
      and gradients.
    - `opf_vlim_hess()` evaluates Hessian of voltage magnitude limit
      constraints.
    - `opf_vref_fcn()` evaluates reference voltage angle equality
      constraints and gradients.
    - `opf_vref_hess()` evaluates Hessian of reference voltage angle
      equality constraints.
    - `opt_model/add_lin_constraint()` to add linear constraints to an
      optimization model.
    - `opt_model/add_nln_constraint()` to add nonlinear constraints to
      an optimization model.
    - `opt_model/init_indexed_name()` to initialize the indices for an
      indexed name set of constraints, costs or variables.
    - `save2psse()` to export a MATPOWER case to PSS/E RAW data format.
    - `savechgtab()` to save change tables, such as those used by
     `apply_changes`, to a file.

#### New Case Files:
- Seven new purely synthetic cases from the ACTIVSg team (**A**SU,
  **C**ornell, **T**exas A&M, U of **I**llinois, and **V**CU -
  **S**ynthetic **g**rids), resulting from work supported by the ARPA-E
  GRID DATA program. *Thanks to Adam Birchfield and the ACTIVSg team.*
  - `case_ACTIVSg200` -- 200-bus Illinois synthetic model
  - `case_ACTIVSg500` -- 500-bus South Carolina synthetic model
  - `case_ACTIVSg2000` -- 2000-bus Texas synthetic model
  - `case_ACTIVSg10k` -- 10,000-bus US WECC synthetic model
  - `case_ACTIVSg25k` -- 25,000-bus US Northeast/Mid-Atlantic synthetic
    model
  - `case_ACTIVSg70k` -- 70,000-bus Eastern US synthetic model
  - `case_SyntheticUSA` -- 82,000-bus continental USA synthetic model
    (aggregation of `case_ACTIVSg70k`, `case_ACTIVSg10k`, and
    `case_ACTIVSg2000`, connected by 9 DC lines)

  Some of these cases also include contingency tables and/or hourly
  load scenarios for 1 year.
  - `contab_ACTIVSg200`
  - `contab_ACTIVSg500`
  - `contab_ACTIVSg2000`
  - `contab_ACTIVSg10k`
  - `scenarios_ACTIVSg200`
  - `scenarios_ACTIVSg2000`
- New RTS-GMLC case from
  [https://github.com/GridMod/RTS-GMLC][3].
  - `case4_RTS_GMLC`
- Six new radial distribution system cases. *Thanks to Mirko Todorovski.*
  - `case4_dist`
  - `case18`
  - `case22`
  - `case69`
  - `case85`
  - `case141`

#### New Documentation:
- Two new Tech Notes, available from [MATPOWER home page][4].
  - B. Sereeter and R. D. Zimmerman, "Addendum to AC Power Flows and
    their Derivatives using Complex Matrix Notation: Nodal Current
    Balance", [*MATPOWER Technical Note 3*][5], April 2018.
  - B. Sereeter and R. D. Zimmerman, "AC Power Flows and their
    Derivatives using Complex Matrix Notation and Cartesian
    Coordinate Voltages", [*MATPOWER Technical Note 4*][6], April 2018.
- LaTeX source code for [MATPOWER User's Manual][7] included in
  `docs/src`, for [MIPS User's Manual][8] in `mips/docs/src` and for
  [MOST User's Manual][9] in `most/docs/src`.

#### Other Improvements:
- Update versions of included packages:
  - MIPS 1.3.
  - MOST 1.0.1.
  - MP-Test 7.0b1.
- Continuous integration testing via GitHub and Travis-CI integration.
- Support added in core optimization model `opt_model` for:
  - general nonlinear constraints
  - general nonlinear costs
  - quadratic costs
- Refactor OPF code to take advantage of new `opt_model` capabilities
  for nonlinear constraints and quadratic and nonlinear costs.
- Derivative functions now support cartesian coordinates for voltage in
  addition to polar coordinates.
- In the Newton power flow, for larger systems use explicit LU
  decomposition with AMD reordering and the 3 output argument form of
  `lu` (to select the Gilbert-Peierls algorithm), resulting in up to a
  2x speedup in MATLAB, 1.1x in Octave. *Thanks to Jose Luis Marín.*
- Support plotting of multiple nose curves in CPF by allowing option
  `cpf.plot.bus` to take on vector values.
- Add line for curtailed load to `case_info()` output.
- Change default implementation of active power line flow constraints
  (`opf.flow_lim = 'P'`) to use flow directly, rather than square of
  flow, which is now a separate option, namely `opf.flow_lim = '2'`.
  *Thanks to Nico Meyer-Huebner*.
- Add `genfuels` and `gentypes` to establish standard set of values for
  optional `mpc.genfuel` and `mpc.gentype` fields for generator fuel
  type and generator unit type, respectively.
- Add support for `gentype` and `genfuel` fields of MATPOWER case
  struct in `extract_islands`, `ext2int`, `int2ext`, `load2disp`  and
  `savecase`.
- Add support for `bus_name` field of MATPOWER case struct to
  `extract_islands`, `ext2int` and `int2ext`.
- Deprecated functions:
    - `d2AIbr_dV2()` -- use `dA2br_dV2()` instead.
    - `d2ASbr_dV2()` -- use `dA2br_dV2()` instead.
    - `opt_model/add_constraints()` -- use the corresponding one of the
      following methods instead: `add_lin_constraint()`,
      `add_nln_constraint()`, or `init_indexed_name()`.
    - `opt_model/add_costs()` -- use the corresponding one of the
      following methods instead: `add_quad_cost()`,  `add_nln_cost()`,
      `add_legacy_cost()`, or `init_indexed_name()`.
    - `opt_model/linear_constraints()` -- use
      `opt_model/params_lin_constraint()` instead.
    - `opt_model/build_cost_params()` -- no longer needed, incorporated
      into `opt_model/params_legacy_cost()`.
    - `opt_model/get_cost_params()` -- use
      `opt_model/params_legacy_cost()` instead.

#### Bugs Fixed:
- Fix bug in conversion of older versions of MATPOWER options.
- Fix bug #4 where some Q limits were not being respected by CPF when
  buses were converted to PQ by initial power flow run.
  *Thanks to Shruti Rao.*
- Fix fatal bug #8 when calling `runcpf` with base and target cases
  with identical load and generation. *Thanks to Felix.*
- Fix fatal bug in `get_losses` when computing derivatives of reactive
  branch injections and fix some related tests.
- Fix #11 fatal error encountered when running `test_matpower` with
  `SDP_PF` and YALMIP installed, but no SDP solver. Now checks for
  availability of SeDuMi, SDP3 or MOSEK before attempting to run
  `SDP_PF` tests that require solving an SDP. *Thanks to Felix.*
- Fix bug #12 where the CPF could terminate early when requesting trace
  of the full curve with P or Q limits enforced, if a limit becomes
  binding at the base case. *Thanks to Felix.*
- Fix bug #13 where setting all buses to type `NONE` (isolated)
  resulted in a fatal error for `ext2int`, `runpf`, `runcpf` and
  `runopf`. *Thanks to SNPerkin.*
- Fix bug #21 where a continuation power flow that failed the first
  corrector step would produce a fatal error. *Thanks to Elis Nycander.*
- Fix bug #23 where the continuation power flow could switch directions
  unexpectedly when the operating point switched from stable to
  unstable manifold or vice-versa after hitting a limit.
  *Thanks to Elis Nycander and Shrirang Abhyankar.*
- Fix bug #26 where, in a continuation power flow, a reactive limit at
  a bus could be detected in error if multiple generators at the bus
  had reactive ranges of very different sizes.
  *Thanks to Elis Nycander and Shrirang Abhyankar.*
- Fix `runpf` handling of case where individual power flow fails during
  Q limit enforcement.

#### Incompatible Changes:
- Move included MATPOWER case files to new `data` subdirectory.
- Turning soft limits on without specifying any parameters explicitly
  in `mpc.softlims` now implements soft limits for all constraints,
  by default, not just branch flow limits. And the format of the input
  parameters in `mpc.softlims` has changed. See `help toggle_softlims`
  or Tables 7-9, 7-10 and 7-11 in the [MATPOWER User's Manual][2] for
  the details.
- Swap the order of the output arguments of `dSbus_dV()`) for polar
  coordinate voltages (angle before magnitude) for consistency.
- Correct signs of phase shifter angles in Polish system cases, since
  they were based on the old sign convention used by MATPOWER prior to
  v3.2 (see change on 6/21/07). Affects the following cases:
  - `case2383wp`
  - `case2736sp`
  - `case2737sop`
  - `case2746wop`
  - `case2746wp`
  - `case3375wp`

  *Thanks to Mikhail Khokhlov and Dr. Artjoms Obusevs for reporting.*
- Remove `nln.mu.l.<name>` and `nln.mu.u.<name>` fields from OPF
  `results` struct. Use `nle.lambda.<name>` and `nli.mu.<name>` fields
  instead for nonlinear constraint multipliers.
- Modify order of default output arguments of `opt_model/get_idx()`.
- Add `mpopt` to input args for OPF `'ext2int'`, `'formulation'`, and
  `'int2ext'` callbacks.


[1]: https://github.com/MATPOWER/matpower/blob/master/CHANGES.md
[2]: https://github.com/MATPOWER/matpower/blob/master/docs/MATPOWER-manual.pdf
[3]: https://github.com/GridMod/RTS-GMLC
[4]: http://www.pserc.cornell.edu/matpower/
[5]: http://www.pserc.cornell.edu/matpower/TN3-More-OPF-Derivatives.pdf
[6]: http://www.pserc.cornell.edu/matpower/TN4-OPF-Derivatives-Cartesian.pdf
[7]: http://www.pserc.cornell.edu/matpower/docs/MATPOWER-manual-7.0b1.pdf
[8]: http://www.pserc.cornell.edu/matpower/docs/MIPS-manual-1.3.pdf
[9]: http://www.pserc.cornell.edu/matpower/docs/MOST-manual-1.0.1.pdf
