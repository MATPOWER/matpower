What's New in MATPOWER 6.0
--------------------------

#### Released Dec 16, 2018

Below is a summary of the changes since version 5.1 of MATPOWER. See the
[`CHANGES.md`][1] file for all the gory details. For release notes for
previous versions, see Appendix H of the [MATPOWER User's Manual][2].

#### New Open Development Model
  - MATPOWER development has moved to GitHub! The [code repository][4] is
    now publicly available to clone and submit pull requests.
  - [Public issue tracker][4] for reporting bugs, submitting patches, etc.
  - Separate repositories for [MOST][5], [MIPS][6], [MP-Test][7],
    all available from the [MATPOWER Development][8] page (a GitHub
    organization).
  - New developer e-mail list ([MATPOWER-DEV-L][9]) to facilitate
    communication between those collaborating on MATPOWER development.

#### New Case Files:
  - Added 9 new case files, 8 cases ranging from 1888 to 6515 buses
    representing the French system, and a 13,659-bus case representing
    parts of the of the European high voltage transmission network,
    stemming from the Pan European Grid Advanced Simulation and State
    Estimation (PEGASE) project. *Thanks again to Cedric Josz and
    colleagues from the French Transmission System Operator.*
  - Added `case145.m`, IEEE 145 bus, 50 generator dynamic test case
    from the [UW Power Systems Test Case Archive]
    (http://www.ee.washington.edu/research/pstca/dyn50/pg_tcadd50.htm).
  - Added `case33bw.m`, a 33-bus radial distribution system from Baran
    and Wu.

#### New Features:
  - [MATPOWER Optimal Scheduling Tool (MOST)][5] is a major new feature,
    implementing a full range of optimal power scheduling problems, from a
    simple as a deterministic, single period economic dispatch problem
    with no transmission constraints to as complex as a stochastic,
    security-constrained, combined unit-commitment and multiperiod OPF
    problem with locational contingency and load-following reserves,
    ramping costs and constraints, deferrable demands, lossy storage
    resources and uncertain renewable generation.
    See [`docs/MOST-manual.pdf`][10] for details.
  - General mechanism for applying modifications to an existing MATPOWER
    case. See `apply_changes()` and `idx_ct()`.
  - Redesigned CPF callback mechanism to handle CPF events such as
    generator limits, nose point detection, etc. Included event log
    in CPF results.
  - Added options `cpf.enforce_p_lims` and `cpf.enforce_q_lims` to
    enforce generator active and reactive power limts in the
    continuation power flow.
  - Added OPF option `opf.use_vg` to provide a convenient way to have
    the OPF respect the generator voltage setpoints specified in the
    gen matrix.
  - Experimental foundation for handling of ZIP load models in power flow
    (Newton, fast-decoupled only), continuation power flow, and optimal
    power flow (MIPS, fmincon, Knitro, IPOPT solvers only). Currently,
    ZIP loads can only be specified on a system-wide basis using the
    experimental options `exp.sys_wide_zip_loads.pw` and
    `exp.sys_wide_zip_loads.qw`.
 - Support for `quadprog()` under GNU Octave.
 - New contributed extras:
    - Plot electrically meaningful drawings of a MATPOWER case using
      `plot_mpc()` in `extras/misc`, *contributed by Paul Cuffe*.
    - Find the maximum loadability limit of a system via an optimal power
      flow and dispatchable loads, using `maxloadlim()` in `extras/maxloadlim`,
      *contributed by Camille Hamon*.
    - Create a quadratically-constrained quadratic programming (QCQP)
      representation of the AC optimal power flow problem using using
      `qcqp_opf()` in `extras/misc`, *contributed by Cedric Josz and
      colleagues*.
  - New functions:
    - `apply_changes()` and `idx_ct()` provide a general mechanism for
      applying modifications to an existing MATPOWER case.
    - `feval_w_path()` evaluates a function located at a specified path,
      outside of the MATLAB path.
    - `mpopt2qpopt()` provides a common interface for creating options
      struct for `mi/qps_matpower()` from a MATPOWER options struct.
  - New function options:
    - Option to call `makeB()`, `makeBdc()`, `makePTDF()`, `scale_load()`,
      and `total_load()` with full case struct (`mpc`) instead of
      individual data matrices (`bus`, `branch`, etc.).
    - `total_load()`, which now computes voltage-dependent load values,
      accepts the values `bus` and `area` as valid values for `load_zone`
      argument.

#### Other Improvements:
  - Changed default solver order for LP, QP, MILP, MIQP problems to move
    Gurobi before CPLEX and BPMPD after OT and GLPK.
  - Added some caching to `mpoption()` and made minor changes to
    `nested_struct_copy()` to greatly decrease the overhead added by
    `mpoption()` when running many small problems.
  - Added option `cpf.adapt_step_damping` to control oscillations in
    adaptive step size control for continuation power flow.
  - Added CPF user options for setting tolerances for target lambda
    detection and nose point detection, `cpf.target_lam_tol` and
    `cpf.nose_tol`, respectively.
  - Added support for MATLAB Optimization Toolbox 7.5 (R2016b).
  - Added support for MOSEK v8.x.
  - Added tests for power flow with `pf.enforce_q_lims` option.
  - Updated network reduction code to handle cases with radially
    connected external buses.
  - Updated versions of `qcqp_opf()` and `qcqp_opf()` in `extras/misc`,
    *from Cedric Josz*.
  - Added "Release History" section to Appendix of manual.
  - Many new tests.

#### Bugs Fixed:
  - Fixed bug in `toggle_dclines()` that resulted in fatal error when used
    with OPF with reactive power costs. *Thanks to Irina Boiarchuk.*
  - Fixed fatal bug in `update_mupq()` affecting cases where `QMIN` is greater
    than or equal to `QC1MIN` and `QC2MIN` (or `QMAX` is less than or equal to
    `QC1MAX` and `QC2MAX`) for all generators. *Thanks Jose Miguel.*
  - Copying a field containing a struct to a non-struct field with
    `nested_struct_copy()` now overwrites rather than causing a fatal error.
  - Fixed a bug in `psse_convert_xfmr()` where conversion of data for
    transformers with CZ=3 was done incorrectly. *Thanks to Jose Marin
    and Yujia Zhu.*
  - Fixed a fatal bug in `psse_convert_xfmr()` affecting transformers with
    CW and/or CZ equal to 1. *Thanks to Matthias Resch.*
  - Fixed a crash in `have_fcn()` caused by changes in OPTI Toolbox v2.15
    (or possibly v2.12)
  - Commented out isolated bus 10287 in `case3375wp.m`.
  - Added code to DC OPF to return `success` = 0 for cases where the matrix
    is singular (e.g. islanded system without slack).
  - Fixed problem in `have_fcn()` where SeDuMi was turning off and leaving
    off all warnings.
  - Fixed shadow prices on variable bounds for AC OPF for fmincon,
    IPOPT, and Knitro.
  - In `savecase()` single quotes are now escaped properly in bus names.
  - Generator capability curve parameters that define a zero-reactive
    power line no longer cause a fatal error.
  - Bad bus numbers no longer cause a fatal error (after reporting the
    bad bus numbers) in `case_info()`.

#### Incompatible Changes:
  - Removed `fairmax()` from the public interface by moving it inside `uopf()`,
    the only place it was used.
  - Removed `cpf.user_callback_args` option and modified
    `cpf.user_callback`.
  - Changed name of `cpf.error_tol` option to `cpf.adapt_step_tol`.


[1]: https://github.com/MATPOWER/matpower/blob/master/CHANGES.md
[2]: https://github.com/MATPOWER/matpower/blob/master/docs/MATPOWER-manual.pdf
[3]: https://github.com/MATPOWER/matpower
[4]: https://github.com/MATPOWER/matpower/issues
[5]: https://github.com/MATPOWER/most
[6]: https://github.com/MATPOWER/mips
[7]: https://github.com/MATPOWER/mptest
[8]: https://github.com/MATPOWER/
[9]: http://www.pserc.cornell.edu/matpower/mailinglists.html#devlist
[10]: https://github.com/MATPOWER/most/blob/master/docs/MOST-manual.pdf
