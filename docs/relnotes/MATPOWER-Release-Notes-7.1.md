What's New in MATPOWER 7.1
--------------------------

#### Released Oct 8, 2020

Below is a summary of the changes since version 7.0 of MATPOWER. See the
[`CHANGES.md`][1] file for all the gory details. For release notes for
previous versions, see Appendix H of the [MATPOWER User's Manual][2].

#### New Features:
- [MP-Opt-Model][3] 3.0 brings many new enhancements to the `opt_model`
  class and the various optimization solver interfaces, including:
  - New unified interface `nlps_master()` for nonlinear programming problems,
    including implementations for [MIPS][4], `fmincon`, IPOPT and
    Artelys Knitro.
  - New unified interface `nleqs_master()` for solving nonlinear equations,
    including implementations for `fsolve`, Newton's method, fast-decoupled
    Newton and Gauss-Seidel.
  - Automatic determination of explicit model type from characteristics of
    variables, costs and constraints.
  - New `solve()` method for `opt_model` to call appropriate solver for the
    respective problem type, including for linear and nonlinear equations.
  - Performance improvements.
  - Support for new solvers ([OSQP][6], `fsolve`, Newton, fast-decoupled
    Newton, Gauss-Seidel) and new versions of CPLEX and Artelys Knitro.
  - Support for modifying parameters of an existing model.
  - Support for extracting specific variables, costs, constraint values
    and shadow prices, etc. from a solved model.

  For more details on improvements related to [MP-Opt-Model][3], see the
  [release notes][7] *(also in Appendix C in the [MP-Opt-Model User's
  Manual][8])* for MP-Opt-Model 0.8, [1.0][9], [2.0][10], [2.1][11], and
  [3.0][12].
- [MP-Test][5] 7.1, with new modular, extensible `have_feature()` function
  for detecting optional functionality. For more details, see the
  [MP-Test 7.1 release notes][13].
- [MIPS][4] 1.4. For details, see the [MIPS 1.4 release notes][14]
  *(also in Appendix C in the [MIPS User's Manual][15])*.
- Support for [OSQP][6] to solve LP and QP problems. Set option
  `opf.dc.solver` to `'OSQP'` to use OSQP to solve the DC OPF. Requires the
  MATLAB interface to OSQP, available from [https://osqp.org][6]. See
  `help mpoption` for more OSQP options.
- Option for `makePTDF()` to compute shift factors very efficiently for
  specific transfers or specific buses only.
- New options:
  - `pf.alg` now accepts `'NR-SP'` as a shortcut for Newton's method with
    the default power/polar formulation, and `'NR-SH'` and `'NR-IH'` for
    those with the hybrid voltage updates.
  - `pf.v_cartesian` now accepts `2` as a valid value to select the Newton
    power flow with hybrid voltage update.
  - `opf.dc.solver` now accepts `'OSQP'` as a valid value if [OSQP][6] is
    installed, to select [OSQP][6] to solve the DC OPF.
  - `osqp.opts` overrides default [OSQP][6] options.
- New functions/methods:
  - MP-Test 7.1
    - `have_feature()` detects availability and version information for
      optional functionality. This is a modular, extensible replacement for
      `have_fcn()` where the detection of a feature named `<tag>` is
      implemented by the function `have_feature_<tag>()`.
  - MP-Opt-Model 3.0
    - `mpopt2nleqopt()` creates or modifies an options struct for
      `nleqs_master()` from a MATPOWER options struct.
    - `mpopt2nlpopt()` creates or modifies an options struct for
      `nlps_master()` from a MATPOWER options struct.
    - `nleqs_master()` provides a single wrapper function for calling any
      of MP-Opt-Model's nonlinear equation solvers.
    - `nlps_master()` provides a single wrapper function for calling any
      of MP-Opt-Model's nonlinear programming solvers.
    - `opt_model.eval_lin_constraint()` method computes the constraint values
      for the full set or an individual named subset of linear constraints.
    - `opt_model.get_soln()` method extracts solved results for a given named
      set of variables, constraints or costs.
    - `opt_model.params_nln_constraint()` method returns parameters for a
      named (and optionally indexed) set of nonlinear constraints.
    - `opt_model.params_nln_cost()` method returns parameters for a named
      (and optionally indexed) set of general nonlinear costs.
    - `opt_model.parse_soln()` method returns a complete set of solution
      vector and shadow price values for a solved model.
    - `opt_model.problem_type()` determines the type of problem based on the
      characteristics of the variables, costs and constraints in the model.
    - `opt_model.set_params()` method modifies parameters for a given named
      set of existing variables, costs, or constraints of an MP-Opt-Model object.
    - `opt_model.solve()` solves the model using `mplinsolve()`,
      `qps_master()`, `miqps_master()`, `nleqs_master()` or `nlps_master()`,
      depending on the problem type.
    - `osqp_options()` initializes options for [OSQP][6] solver.
    - `osqpver()` returns/displays version information for [OSQP][6].
    - ... plus `nleqs_core()`, `nleqs_fd_newton()`, `nleqs_fsolve()`,
      `nleqs_gauss_seidel()`, `nleqs_newton()`, `nlps_fmincon()`,
      `nlps_ipopt()`, `nlps_knitro()`, `opt_model.copy()`,
      `opt_model.is_mixed_integer()`, and `qps_osqp()`.
- `nlpopf_solver()` implements the AC OPF solver based on `opt_model`'s new
  `solve()` method, and replaces the individual solver-specific functions
  `fmincopf_solver()`, `ipoptopf_solver()`, `ktropf_solver()`, and 
  `mipsopf_solver()`.
- ... plus 45 individual feature detection functions for `have_feature()`
  *(10 in MATPOWER proper, plus 2 in MP-Test, 4 in MIPS, and 29 in MP-Opt-Model
  (see Table A-7 in the [MP-Opt-Model User's Manual][8] for details).*.

#### New Case Files:
- Nineteen new distribution system cases \[[ref 1][16], [ref 2][17]\].
  *Thanks to Houssem Bouchekara, et. al.*
  - `case10ba` -- 10-bus radial system from Baghzouz and Ertem
  - `case12da` -- 12-bus radial system from Das, Nagi, and Kothari
  - `case15da` -- 15-bus radial system from Das, Kothari, and Kalam
  - `case15nbr` -- 15-bus radial system from Battu, Abhyankar, Senroy
  - `case16am` -- 16-bus radial system from Das, Kothari, and Kalam
  - `case16ci` -- 16-bus system from Civanlar, Grainger, Yin, and Lee
  - `case17me` -- 17-bus radial system from Mendoza, Morales, Lopez, et. al.
  - `case18nbr` -- 18-bus radial system from Battu, Abhyankar, Senroy
  - `case28da` -- 28-bus radial system from Das, Nagi, and Kothari
  - `case33mg` -- 33-bus radial system from Kashem, et. al.
  - `case34sa` -- 34-bus radial system from Salama and Chikhani
  - `case38si` -- 38-bus radial system from Singh and Misra
  - `case51ga` -- 51-bus radial system from Gampa and Das
  - `case51he` -- 51-bus radial system from Hengsritawat, et. al.
  - `case70da` -- 70-bus system from Das
  - `case74ds` -- 74-bus radial system from Myint and Naing
  - `case94pi` -- 94-bus radial system from Pires, Antunes and Martins
  - `case118zh` -- 118-bus radial system from Zhang, Fu and Zhang
  - `case136ma` -- 136-bus radial system from Mantovani, Casari and Romero

#### New Documentation:
- New [MP-Opt-Model User's Manual][8], included in `mp-opt-model/docs`.

#### Other Improvements:
- Refactor `opt_model` class to inherit from new abstract base class
  `mp_idx_manager` which can be used to manage the indexing of other sets
  of parameters, etc. in other contexts.
- Add to `opt_model.eval_nln_constraint()` method the ability to compute
  constraints for a single named set.
- Significant performance improvement for some problems when constructing
  sparse matrices for linear constraints or quadratic costs (e.g. during
  problem setup in MOST). *Thanks to Daniel Muldrew.*
- Move original implementation of Newton power flow for cartesian
  voltage representations to `newtonpf_S_hybrid()` and `newtonpf_I_hybrid()`,
  accessible by setting the `pf.v_cartesian` option to 2. The `pf.alg` option
  can also be set to `'NR-SH'` or `'NR-IH'`, respectively, as shortcuts to
  select these formulations.
- Improve robustness of these hybrid Newton power flow formulations to avoid
  ill-conditioning when the real part of the voltage approaches zero.
- Redesign implementation of Newton power flow for cartesian voltage
  representations to use standard Newton method with explicit voltage
  magnitude constraints at PV buses. This method is slower (larger number of
  equations and unknowns) but appears to be more reliable in some cases than
  the original implementation, which was a hybrid approach by Sereeter that
  used a modified Newton method to compute polar voltage updates using a
  modified cartesian Jacobian.
- Modify voltage limit constraints for cartesian AC OPF formulation to use
  voltage squared, resulting in simpler derivatives.
- Performance improvement for `makePTDF()` for large cases, e.g. more
  than 70% improvement for `case9241pegase`.
- Significant performance improvement for CPLEX on small problems by
  eliminating call to `cplexoptimset()`, which was a huge bottleneck.
- Convert `dcopf_solver()` to use the new `solver()` method of `opt_model`
  instead of calling `qps_matpower()` manually.
- Convert `opf_execute()` to use `nlpopf_solver()`, based on the new
  `solver()` method of `opt_model`, for AC OPF when using `fmincon`, IPOPT,
  Artelys Knitro, or [MIPS][4].
- Reduce memory usage in `modcost()` for very large systems.
  *Thanks to Christoph Funke.*
- Deprecated functions:
  - `have_fcn()` -- use `have_feature()` from [MP-Test][5] instead.
  - `qps_matpower()` -- use `qps_master()` from [MP-Opt-Model][3] instead.
  - `miqps_matpower()` -- use `miqps_master()` from [MP-Opt-Model][3] instead.
- Removed functions:
  - `fmincopf_solver()` -- functionality now covered by `nlpopf_solver()`.
  - `ipoptopf_solver()` -- functionality now covered by `nlpopf_solver()`.
  - `ktropf_solver()` -- functionality now covered by `nlpopf_solver()`.
  - `mipsopf_solver()` -- functionality now covered by `nlpopf_solver()`.

#### Bugs Fixed:
- For `opt_model`, incorrect evaluation of constant term has been fixed for
  vector valued quadratic costs with constant term supplied as a vector.
- Calling `opt_model.params_var()` method with empty `idx` no longer results
  in fatal error.
- Fix bug in `scale_load()` where only one of multiple dispatchable loads at
  a bus would be scaled.
  *Thanks to Christoph Funke.*
- Fix bug #77 where incorrect indexing could cause fatal error in OPF with
  additional nonlinear constraints.
  *Thanks to Sergio Garcia.*
- Fix OPF issue #71 for IPOPT and Artelys Knitro where structure of Jacobian
  and/or Hessian could change from the structure provided (i.e. elements with
  value of zero were being eliminated from the sparsity structure).
  *Thanks to Drosos Kourounis.*
- Artelys Knitro 12.1 compatibility fix.
- Fix CPLEX 12.10 compatibility issue #90.
- Fix bug #89 where running a power flow with `pf.enforce_q_lims` enabled and
  voltage dependent ZIP loads produced incorrect results.
  *Thanks to Florian.*
- Fix issue with missing objective function value from `miqps_mosek()` and
  `qps_mosek()` when return status is "Stalled at or near optimal solution."
- Fix bug orginally in `ktropf_solver()` (code now moved to `nlps_knitro()`)
  where Artelys Knitro was still using `fmincon` options.

#### Incompatible Changes:
- Modify order of default output arguments of `opt_model.get_idx()` (again),
  removing the one related to legacy costs.
- MP-Opt-Model has renamed the following functions and modified the order of
  their input args so that the MP-Opt-Model object appears first. Ideally,
  these would be defined as methods of the `opt_model` class, but Octave 4.2
  and earlier is not able to find them via a function handle (as used in the
  `solve()` method) if they are inherited by a subclass.
  - `opf_consfcn()` ––> `nlp_consfcn()`
  - `opf_costfcn()` ––> `nlp_costfcn()`
  - `opf_hessfcn()` ––> `nlp_hessfcn()`
- Update `case18`, `case22`, `case69`, `case85` and `case141` to more closely
  match data from original papers, thanks in part to case files contributed
  by Houssem Bouchekara, et al. Solutions for updated cases may not match
  exactly. See help text in case files for details.


[1]: ../../CHANGES.md
[2]: ../MATPOWER-manual.pdf
[3]: https://github.com/MATPOWER/mp-opt-model
[4]: https://github.com/MATPOWER/mips
[5]: https://github.com/MATPOWER/mptest
[6]: https://osqp.org
[7]: https://github.com/MATPOWER/mp-opt-model/blob/master/docs/relnotes/
[8]: https://matpower.org/docs/MP-Opt-Model-manual-3.0.pdf
[9]: https://github.com/MATPOWER/mp-opt-model/blob/master/docs/relnotes/MP-Opt-Model-Release-Notes-1.0.md
[10]: https://github.com/MATPOWER/mp-opt-model/blob/master/docs/relnotes/MP-Opt-Model-Release-Notes-2.0.md
[11]: https://github.com/MATPOWER/mp-opt-model/blob/master/docs/relnotes/MP-Opt-Model-Release-Notes-2.1.md
[12]: https://github.com/MATPOWER/mp-opt-model/blob/master/docs/relnotes/MP-Opt-Model-Release-Notes-3.0.md
[13]: https://github.com/MATPOWER/mptest/blob/master/docs/relnotes/MP-Test-Release-Notes-7.1.md
[14]: https://github.com/MATPOWER/mips/blob/master/docs/relnotes/MIPS-Release-Notes-1.4.md
[15]: https://matpower.org/docs/MIPS-manual-1.4.pdf
[16]: https://doi.org/10.18799/24056537/2019/3/196
[17]: https://doi.org/10.36227/techrxiv.12578648.v1
