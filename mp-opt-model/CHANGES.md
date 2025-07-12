Change history for MP-Opt-Model
===============================


Version 5.0 - *Jul 12, 2025*
----------------------------

#### 7/12/25
  - Release 5.0.
  - Move examples from `lib/t` to their own `examples` directory.

#### 7/3/25
  - Add support for Artelys Knitro 15.x which required changes to the
    prior options handling.
  - **INCOMPATIBLE CHANGE:** The `knitro_opts` field of the `opt` input
    to `nlps_master()` and `nlps_knitro()` and the `solve()` method of
    `opt_model` has been redesigned. It is now a raw Artelys Knitro options
    struct, so the `opts`, `tol_x` and `tol_f` fields are no longer valid.
    For `tol_x` and `tol_f`, use `xtol` and `ftol`, and the contents of
    `opts` should be placed directly in the top level of the `knitro_opts`
    field.
  - **INCOMPATIBLE CHANGE:** Remove support for older versions of Knitro,
    including all references to `ktrlink` for pre-v9 versions. Currently
    supports Artelys Knitro version 13.1 and later.

#### 6/21/25
  - Add new `mp.struct2object()` function to convert a struct back to the
    object from which it was created. Helps with workarounds to the fact
    that Octave still (as of 10.x) does not support saving/loading of
    classdef objects. This function allows objects to implement
    `to_struct()` and `from_struct()` methods to facilitate the process.
  - Add `to_struct()` and `from_struct()` methods to `mp.opt_model`,
    `mp.set_manager`, and legacy `opt_model` classes, to facilitate
    trivial conversion between objects, which Octave cannot save/load,
    and structs, which it can.

#### 6/10/25
  - Fix handling of scalar inputs for vector parameters when adding an
    empty set of variables or linear/quadratic constraints. Now
    properly "expands" them to an empty vector ([issue #16][18]).
  - Update handling by `mp.sm_quad_cost` of constant term for
    quadratic costs ([issue #15][17]):
    - When `H` is empty, a scalar `k` will no longer be expanded
      _(implicitly)_ to a vector, rather it will result in a scalar
      cost set.
    - When `H` is a vector and `k` is a scalar, `k` will be expanded
      _explicitly_ to a vector.
    The `mp.sm_quad_cost_legacy` class is unchanged, so this change
    affects only `mp.opt_model`, not the legacy `opt_model`.

#### 5/23/25
  - Add support to `qps_master()` and `miqps_master()` for the
    open-source [HiGHS][15] solver for LP, QP and MILP problems,
    including functions `miqps_highs()`, `qps_highs()`,
    `have_feature_highs()`, `highsver()`, and `highs_options()`. Uses
    the [HiGHSMEX][15] interface for HiGHS.
  - Fix that using `miqps_master()` with `'DEFAULT'` solver to solve an
    LP/QP problem without a MILP/MIQP solver available incorrectly
    threw a fatal error stating there was no solver available.

#### 5/7/25
  - Fix crash in `mp.sm_quad_cost.display_soln()` when linear
    coefficient is empty.

#### 4/30/25
  - Add support for quadratic constraint and QCQP (quadratically-
    constrained quadratic programming) problems.
    *Thanks to Wilson Gonzalez Vanegas.*
    - Add new top-level wrapper function `qcqps_master()` to provide
      a standard unified interface quadratically-constrained quadratic
      programming (QCQP) problems, with the ability to provide solver-
      specific input options.
    - Add `qcqps_gurobi()`, `qcqps_knitro()`, `qcqps_nlps()` with
      interface that matches `qcqps_master()` to handle implementation
      for Gurobi, Artelys Knitro, `fmincon`, IPOPT, and MIPS  solvers.
    - Add `mpopt2qcqpopt()` to set up an options struct for
      `qcqps_master()` based on a MATPOWER options struct.
    - Add `mp.sm_quad_constraint` class to handle quadratic
      constraints in `opt_model`.
    - Add `qcn` property to `opt_model` for quadratic constraints
      and automatic detection of a new `QCQP` problem type that is
      sent to `qcqps_master()` to solve.

#### 3/19/25
  - Update `miqps_<solver>()` functions to avoid changing MIP solution
    values in price computation stage. It was rounding integer variables,
    potentionally causing a small discrepancy between the objective
    value reported by the solver and the value obtained by computing
    directly from the returned solution x.

#### 3/8/25
  - Update `mosek_options()` for MOSEK 11.x compatibility.

#### 3/1/25
  - Add support to `qps_master()` for [Artelys Knitro][14] solver for
    LP and QP problems, including functions `qps_knitro()`,
    `knitrover()`, and `artelys_knitro_options()`.
    *Thanks to Wilson Gonzalez Vanegas.*

#### 2/28/25
  - Throw error if `opt_model.parse_soln()` is called for an unsolved
    model.

#### 10/29/24
  - Add `relax_integer` option for `opt_model.solve()`. Set to true to easily
    solve the LP or QP relaxation of a mixed integer LP or QP.

#### 8/17/24
  - Add `mp.set_manager_opt_model` base class to handle common `opt_model`
    functionality, such as handling solutions, for the individual field
    object subclasses.
  - Refactor `opt_model` to move lots of functionality into new
    `mp.set_manager_opt_model` subclasses:
    - `mp.sm_lin_constraint` - set manager class for linear constraints
    - `mp.sm_nln_constraint` - set manager class for nonlinear constraints
    - `mp.sm_nln_cost` - set manager class for general nonlinear costs
    - `mp.sm_quad_cost` - set manager class for quadratic costs
    - `mp.sm_variable` - set manager class for variables
  - Deprecate the following `opt_model` methods in favor of methods of the
    individual `mp.set_manager` objects contained by the `opt_model` object:
    - `add_named_set()` --> use `mp.set_manager.add()`
    - `describe_idx()` --> use `mp.set_manager.describe_idx()`
    - `getN()` --> use `mp.set_manager.get_N()`
    - `init_indexed_name()` --> use `mp.set_manager.init_indexed_name()`
    - `set_type_idx_map()` --> use `mp.set_manager.set_type_idx_map()`
    - `add_lin_constraint()` --> use `mp.sm_lin_constraint.add()`
    - `add_nln_constraint()` --> use `mp.sm_nln_constraint.add()`
    - `add_nln_cost()` --> use `mp.sm_nln_cost.add()`
    - `add_quad_cost()` --> use `mp.sm_quad_cost.add()`
    - `add_var()` --> use `mp.sm_variable.add()`
    - `eval_lin_constraint()` --> use `mp.sm_lin_constraint.eval()`
    - `eval_nln_constraint()` --> use `mp.sm_nln_constraint.eval()`
    - `eval_nln_constraint_hess()` --> use `mp.sm_nln_constraint.eval_hess()`
    - `eval_nln_cost()` --> use `mp.sm_nln_cost.eval()`
    - `eval_quad_cost()` --> use `mp.sm_quad_cost.eval()`
    - `init_indexed_name()` --> use `mp.set_manager.init_indexed_name()`
    - `params_lin_constraint()` --> use `mp.sm_lin_constraint.params()`
    - `params_nln_constraint()` --> use `mp.sm_nln_constraint.params()`
    - `params_nln_cost()` --> use `mp.sm_nln_cost.params()`
    - `params_quad_cost()` --> use `mp.sm_quad_cost.params()`
    - `params_var()` --> use `mp.sm_variable.params()`
    - `set_params()` --> use `mp.set_manager.set_params()`
    - `varsets_cell2struct()` --> use `mp.sm_variable.varsets_cell2struct()`
    - `varsets_idx()` --> use `mp.sm_variable.varsets_idx()`
    - `varsets_len()` --> use `mp.sm_variable.varsets_len()`
    - `varsets_x()` --> use `mp.sm_variable.varsets_x()`
  - **INCOMPATIBLE CHANGE:** Move parsed solution from `soln` property
    of `opt_model` object to `soln` property of individual child
    `mp.set_manager_opt_model` objects. Currently it is still available at
    the original location, but this is now deprecated.

#### 7/10/24
  - Add option for `opt_model.display_soln()` to print to a file.

#### 6/24/24
  - Add `mp.set_manager` class to encapsulate `mp_idx_manager` functionality
    into an individual field object representing a specific set type, rather
    than in the container class.
  - Refactor `mp_idx_manager` to use new `mp.set_manager` class.


Version 4.2 - *May 10, 2024*
----------------------------

#### 5/10/24
  - Release 4.2.

#### 4/23/24
  - Fix bug in test `t_opt_model()` for MATLAB R2011b and earlier.

#### 4/22/24
  - Fix false positive in `have_feature_fsolve()` in case where the file is
    present, but without a valid license.

#### 4/5/24
  - Fix bug in `miqps_mosek()` where the lower and upper bounds of binary
    variables got overwritten with 0 and 1, respectively, effectively
    relaxing any potentially tighter bounds provided as input.
  - Clear cached parameters after updating linear constraints via
    `opt_model.set_params()`.
  - Add caching of aggregate output parameters in `opt_model.params_var()`.

#### 3/26/24
  - Add Sphinx-based [Reference documentation][13].

#### 3/21/24
  - Add to the `parse_soln()` method of `opt_model` an optional `stash`
    input argument that, if present and true, causes the parsed solution
    to be stored back in the object, as the `solve()` method was already
    doing when `opt.parse_soln` is true.
  - Add new method `has_parsed_soln()` to `opt_model` to check for
    availability of a parsed solution in the model.

#### 3/1/24
  - Update `have_feature_ipopt()` to recognize IPOPT MEX installations from
    Enrico Bertolazzi's [mexIPOPT][12], which include MEX files that have
    been renamed to architecture-specific names along with an `ipopt.m`
    wrapper function to call the appropriate one.
    *Thanks to Carlos Murillo-Sánchez.* 
    _**Note:** While MP-Opt-Model no longer requires this, my recommendation
     is still to simply rename the MEX file to `ipopt.<mexext>`, with the
     appropriate architecture-specific extension, and delete the
     unnecessary `ipopt.m` entirely._

#### 1/31/24
  - Add `convert_lin_constraint()` and `convert_lin_constraint_multipliers()`
    functions to eliminate code duplication for common task of converting
    linear constraints and their multipliers between a single set of
    doubly-bounded inequality constraints and separate sets of equality
    and upper-bounded inequality constraints.
  - Change solver for CPLEX price computation stage in `miqps_cplex()` from
    primal simplex to dual simplex (probably no impact, except it was a
    simple way to get a newly failing test in another project to pass again).

#### 12/8/23
  - Always skip price computation stage in `miqps_<solver>()` functions for
    pure (as opposed to mixed) integer problems.

#### 11/29/23
  - Add support to `miqps_master()` for calling `miqps_<my_solver>()` by
    setting `opt.alg` to `'<MY_SOLVER>'` to allow for handling custom
    MILP/MIQP solvers.

#### 11/8/23
  - Add support to `nlps_master()` for calling `nlps_<my_solver>()` by setting
    `opt.alg` to `'<MY_SOLVER>'` to allow for handling custom NLP solvers.

#### 11/6/23
  - Add to `opt_model.add_lin_constraint()` the option to provide/store
    the transpose of the `A` matrix instead of the original. This can
    potentially save significant memory for sparse matrices with many more
    columns than rows. E.g. storage constraints in [MOST][11] for 8760 hour
    planning horizon.

#### 10/13/23
  - Update `opt_model.display_soln()` to avoid displaying an infinite
    average for quadratic costs when corresponding quantity is zero.

#### 9/13/23
  - Clear cached parameters after updating quadratic costs via
    `opt_model.set_params()`.

#### 9/12/23
  - Tweak some MI/QPS solver tests to make them more robust for `'DEFAULT'`
    solver under different environments.

#### 3/27/23
  - Update for compatibility with MATLAB R2023a (Optimization Toolbox 9.5)
    which removed `x0` as a valid input to `linprog`.

#### 3/7/23
  - Add `opt_model.is_solved()` method to determine if the model
    has been solved and contains the solution.
  - Add `opt_model.display_soln()` method to display the results
    of a solved model, including values, bounds and shadow prices for
    variables and linear constraints, values and shadow prices for
    nonlinear constraints, and individual cost components.


Version 4.1 - *Dec 13, 2022*
----------------------------

#### 12/13/22
  - Release 4.1.

#### 11/14/22
  - Allow for Gurobi's `Method` option to be set to 5, for
    _deterministic concurrent simplex_.

#### 10/27/22
  - Add `runtime` field to `output` argument of `qps_glpk()` and
    `qps_mosek()`.
  - Add support to `qps_master()` for calling `qps_<my_solver>()` by setting
    `opt.alg` to `'<MY_SOLVER>'` to allow for handling custom LP/QP solvers.

#### 7/5/22
  - Update for compatibility with Artelys Knitro 13.1 and later.

#### 4/8/22
  - Relax some test tolerances to prevent failure with Gurobi 9.5.

#### 3/3/22
  - Add elapsed time in seconds to results of the `solve()` method of
  `opt_model`, returned in `om.soln.output.et`.

#### 2/14/22
  - Skip some `fmincon` tests using interior point algorithm with finite
    difference Hessian that began failing in MATLAB R2021b.


Version 4.0 - *Oct 18, 2021*
----------------------------

#### 10/18/21
  - Release 4.0.

#### 5/8/21
  - Add several enhancements to `mp_idx_manager/set_type_idx_map()`.
    - Speed improvements
    - Add idxs = [] option to specify "all"
    - Add optional `group_by_name` input argument.

#### 5/3/21
  - Add support for new class of problems - parameterized nonlinear
    equations (PNE). Either create a model with only equality constraints
    (no inequalities or costs) and with number of variables equal to 1 more
    than number of constraints, _or_ call `pnes_master()` directly.
    See Section 4.5 of User's Manual for details.
    - Predictor/corrector numerical continuation method for tracing
      solution curves for PNE problems.
    - Plotting of solution curves.
    - User-defined event functions and callback functions.
    - Warm-start capabilities.

    *Thanks to Shrirang Abhyankar and Alexander Flueck for contributions to
    this feature.*
  - Add functions:
    - `pnes_master()`
    - `pne_callback_default()`
    - `pne_callback_nose()`
    - `pne_callback_target_lam()`
    - `pne_detect_events()`
    - `pne_detected_event()`
    - `pne_event_nose()`
    - `pne_event_target_lam()`
    - `pne_pfcn_arc_length()`
    - `pne_pfcn_natural()`
    - `pne_pfcn_pseudo_arc_length()`
    - `pne_register_callbacks()`
    - `pne_register_events()`
    - `mpopt2pneopt()`

#### 1/21/21
  - Refactor `describe_idx()` into a new method, `set_type_idx_map()`,
    that returns in information in a programmatically usable form, and
    an updated `describe_idx()` that calls the new method, then formats
    the results in the expected char array(s).

#### 1/5/21
  - Calling the `problem_type()` or `is_mixed_integer()` method on an
    empty model no longer cause a fatal error.

#### 12/16/20
  - Update to use labels from `set_types` property as headers for
    `opt_model.display()` to simplify things and facilitate use by
    subclasses.

#### 12/1/20
  - Move `init_set_types()` call out of `opt_model` constructor to avoid
    complexity of other workarounds for [bug in Octave 5.2 and earlier][10]
    related to inheritance of methods called during construction.
    Now `init_set_types()` is called automatically as needed after
    object construction and before object use in `add_var()`, `display()`
    and `init_indexed_name()`.

#### 11/24/20
  - Add optional threshold for detecting failure of LEQ solve, by setting
    the `leq_opt.thresh` option. If the absolute value of any element of
    the solution vector exceeds the threshold, `exitflag` is set to 0,
    indicating failure.

#### 11/17/20
  - Fix examples of `om.set_params()` usage in manual and help text.


Version 3.0 - *Oct 8, 2020*
---------------------------

#### 10/8/20
  - Release 3.0.

#### 9/23/20
  - Move `have_feature()` to [MP-Test][8] and respective feature detection
    functions to [MP-Test][8], [MIPS][9], and [MATPOWER][1].
    - MP-Test
      - `have_feature()`
      - `have_feature_matlab()`
      - `have_feature_octave()`
    - MIPS
      - `have_feature_lu_vec()`
      - `have_feature_pardiso_legacy()`
      - `have_feature_pardiso_object()`
      - `have_feature_pardiso()`
    - MATPOWER
      - `have_feature_e4st()`
      - `have_feature_minopf()`
      - `have_feature_most()`
      - `have_feature_pdipmopf()`
      - `have_feature_regexp_split()`
      - `have_feature_scpdipmopf()`
      - `have_feature_sdp_pf()`
      - `have_feature_smartmarket()`
      - `have_feature_syngrid()`
      - `have_feature_tralmopf()`

#### 9/18/20
  - Add `have_feature()` as a modular, extensible alternative
    to `have_fcn()`, where the detection of a feature named
    `<tag>` is implemented by the function `have_feature_<tag>()`.
  - Make `have_fcn()` a simple wrapper around the new `have_feature()`.

#### 9/16/20
  - Add `set_params()` method to modify parameters of existing
    variables, costs and constraints of an MP-Opt-Model object.
  - Calling `params_var()` method with empty `idx` no longer results
    in fatal error.
  - For `opt_model`, fixed incorrect evaluation of constant term in
    vector valued quadratic costs with constant term supplied as a
    vector.
  - Simplified logic to determine whether a quadratic cost for an
    MP-Opt-Model object is vector vs. scalar valued. If the quadratic
    coefficient is supplied as a matrix, the cost is scalar varied,
    otherwise it is vector valued.

#### 9/14/20
  - Allow `v0`, `vl`, and `vu` inputs to `opt_model.add_var()` method,
    and `l` and `u` inputs to `opt_model.add_lin_constraint()` to
    be scalars that get expanded automatically to the appropriate
    vector dimension.

#### 9/11/20
  - Add `get_soln()` method to `opt_model` for extracting solved
    results for a given named set of variables, constraints or costs.
  - Add `parse_soln()` method which returns a struct with a complete
    set of solution vector and shadow price values for a solved model.

#### 9/10/20
  - Add caching of problem_type() return value.

#### 9/1/20
  - Add support for OSQP solver from [https://osqp.org][7] for LP
    and QP problems, including functions `qps_osqp()`, `osqpver()`,
    and `osqp_options()`.

#### 8/31/20
  - Save the results of `solve()` method to the `soln` field of the
    MP-Opt-Model object.

#### 8/28/20
  - Add `eval_lin_constraint()` method to evaluate the constraint
    values for the full set or an individual named subset of linear
    constraints.

#### 8/27/20
  - Starting point supplied to `solve()` via `opt.x0` is no longer
    ignored for nonlinear equations.


Version 2.1 - *Aug 25, 2020*
----------------------------

#### 8/25/20
  - Release 2.1.
  - Add core NLEQ solver function `nleqs_core()` with arbitrary,
    user-defined update function, used to implement Gauss-Seidel and
    Newton solvers, `nleqs_gauss_seidel()` and `nleqs_newton()`.

#### 8/20/20
  - Add linear equation (`'LEQ'`) problem type for models with equal
    number of variables and linear equality constraints, no costs,
    and no inequality or nonlinear equality constraints. Solved via
    `mplinsolve()`.
  - The `solve()` method of `opt_model` can now automatically handle
    mixed systems of equations, with both linear and nonlinear equality
    constraints.

#### 7/30/20
  - Add fast-decoupled Newton's and Gauss-Seidel solvers for nonlinear
    equations. Use `alg = 'FD'` and `alg = 'GS'` with `nleqs_master()`.
    See also `nleqs_fd_newton()` and `nleqs_gauss_seidel()`.

#### 7/29/20
  - Allow `solve()` method to pass along number of requested output
    arguments `*_master()` solver functions.
  - **INCOMPATIBLE CHANGE:** In `output` return value from
    `nleqs_newton()`, changed the `normF` field of `output.hist` to
    `normf`, for consistency in using lowercase `f` everywhere.


Version 2.0 - *Jul 8, 2020*
---------------------------

#### 7/8/20
  - Release 2.0.

#### 7/3/20
  - Add to `eval_nln_constraint()` method the ability to compute
    constraints for a single named set.

#### 7/2/20
  - Skip evaluation of gradient if `eval_nln_constraint()` is called
    with a single output argument.
  - Add `params_nln_constraint()` method, and add documentation to the
    manual for it and `params_nln_cost()`.

#### 7/1/20
  - Add `mpopt2nleqopt()` to create or modify an `nleqs_master()`
    options struct from a MATPOWER options struct.
  - Add table of valid `have_fcn()` input tags to User's Manual.

#### 6/24/20
  - Add support for nonlinear equations (NLEQ) to `opt_model`. For
    problems with only nonlinear equality constraints and no costs, the
    `problem_type()` method returns `'NLEQ'` and the `solve()` method
    calls `nleqs_master()` to solve the problem.
  - Add tests for solving LP/QP, MILP/MIQP, NLP and NLEQ problems via
    `opt_model.solve()`.

#### 6/17/20
  - Add `nleqs_master()` function as unified interface for solving
    nonlinear equations, including implementations for `fsolve` and
    Newton's method in functions `nleqs_fsolve()` and `nleqs_newton()`,
    respectively.

#### 6/16/20
  - Add new `'fsolve'` tag to `have_fcn()` to check for availability of
    `fsolve()` function.

#### 5/11/20
  - Remove redundant MIPS tests from `test_mp_opt_model`.


Version 1.0 - *May 8, 2020*
---------------------------

#### 5/8/20
  - Release 1.0.

#### 5/7/20
  - Add MP-Opt-Model User's Manual in `docs`, with LaTeX source in
    `docs/src`.
  - Add usage examples to `README.md`.

#### 4/30/20
  - Add `README.md`, `CHANGES.md`, `AUTHORS`, `CONTRIBUTING.md`, `LICENSE`.
  - Refactor `opt_model` class to inherit from new abstract base class
    `mp_idx_manager`, which can be used to manage the indexing of other sets of
    parameters, etc. in other contexts.


Version 0.8 - *Apr 29, 2020*
----------------------------

#### 4/29/20
  - Version 0.8.
  - Add `mpomver()` to define MP-Opt-Model version information.
  - **INCOMPATIBLE CHANGE:** Renamed the following functions and
    modified the order of their input args so that the MP-Opt-Model
    object appears first. Ideally, these would be defined as methods
    of the `opt_model` class, but Octave 4.2 and earlier is not
    able to find them via a function handle (as used in the `solve()`
    method) if they are inherited by a subclass.
    - `opf_consfcn()` --> `nlp_consfcn()`
    - `opf_costfcn()` --> `nlp_costfcn()`
    - `opf_hessfcn()` --> `nlp_hessfcn()`
  - Add GitHub Actions CI workflow and [Travis-CI][3] configuration.
  - Add `test_mp_opt_model()` to run all tests.
  - Remove MATPOWER dependencies.
  - Move code related to solver interfaces, `opt_model` and a
    few other functions like `have_fcn()` and `nested_struct_copy()`
    from main [MATPOWER][1] repository to new [MP-Opt-Model][2]
    repository.


⬆ _In [MP-Opt-Model][2] repository_ ⬆

-----------------------------------

⬇ _In [MATPOWER][1] repository_ ⬇


#### 4/28/20
  - Move deprecated `opt_model` methods and code related to legacy
    user-defined OPF costs from `opt_model` to `opf_model`.
  - **INCOMPATIBLE CHANGE:** Modify order of default output arguments of
    `opt_model.get_idx()` (again), removing the one related to legacy
    costs.

#### 3/18/20
  - Add `nlpopf_solver()` based on the new `solver()` method of
    `opt_model`. This single function replaces `mipsopf_solver()`,
    `fmincopf_solver()`, `ipoptopf_solver()`, and `ktropf_solver()`.
  - Convert `dcopf_solver()` to use the new `solver()` method of
    `opt_model` instead of calling `qps_matpower()` manually.
  - Add new top-level wrapper function `nlps_matpower()` to provide
    a standard interface for MATPOWER's nonlinear program (NLP)
    solvers (`fmincon`, IPOPT, Artelys Knitro, and MIPS), with
    calling syntax similar to `mips()`. It includes the ability to
    pass in solver-specific input options.
  - Add `nlps_fmincon()`, `nlps_ipopt()` and `nlps_knitro()`, with
    interface that matches `nlps_matpower()` to handle implementation
    for `fmincon`, IPOPT, and Artelys Knitro solvers, respectively.
  - Add `mpopt2nlpopt()` to set up an options struct for
    `nlps_matpower()` based on a MATPOWER options struct.
  - Add three new methods to `opt_model` class:
    - `is_mixed_integer()` - returns true if the model includes any binary
      or integer variables
    - `problem_type()` - returns one of the following strings, based on
      the characteristics of the variables, costs and constraints in the
      model:
      - `'NLP'` - nonlinear program
      - `'LP'` - linear program
      - `'QP'` - quadratic program
      - `'MILP'` - mixed-integer linear program
      - `'MIQP'` - mixed-integer quadratic program
    - `solve()` - solves the model using `qps_matpower()`,
      `miqps_matpower()`, or `nlps_matpower()`, depending on the problem
      type (`'MINLP'` problems are not yet implemented)

#### 3/12/20
  - Fix bug in `ktropf_solver()` where Artelys Knitro was still using
    `fmincon` options.

#### 3/6/20
  - Fix issue with missing objective function value from `miqps_mosek()`
    and `qps_mosek()` when return status is "Stalled at or near optimal
    solution."

#### 3/4/20
  - Remove unused input arguments from `opf_consfcn()` and `opf_hessfcn()`.

#### 2/27/20
  - Add `copy()` method to `opt_model` class to get around issues
    with inheritance in constructors that was preventing copy constructor
    from working in Octave 5.2 and earlier (see also [Octave bug
    52614](https://savannah.gnu.org/bugs/?52614).

#### 2/26/20
  - Significant performance improvement for CPLEX on small problems by
    eliminating call to `cplexoptimset()`, which was a huge bottleneck.
  - Fix CPLEX 12.10 compatibility [issue #90][6].

#### 2/18/20
  - Artelys Knitro 12.1 compatibility fix.

#### 8/15/19
  - Improve performance of `opt_model.add_named_set()`.
    (See [issue #79][5].)
    *Thanks to Baraa Mohandes.*
  - Refactor code in `opt_model.params_lin_constraint()` and
    `opt_model.params_quad_cost()` to speed up sparse matrix construction
    when there are lots of constraint or cost sets. Results in significant
    speedups for some problems during problem setup in MOST.
    (See [pull request #70][4].)
    *Thanks to Daniel Muldrew.*


Version 0.7.0 - *Jun 20, 2019*
------------------------------

#### 6/20/19
  - This change history begins with the code that was part of the
    MATPOWER 7.0 release, which is tagged as version 0.7.0 in the
    MP-Opt-Model repository.

----
[1]: https://github.com/MATPOWER/matpower
[2]: https://github.com/MATPOWER/mp-opt-model
[3]: https://travis-ci.org
[4]: https://github.com/MATPOWER/matpower/pull/70
[5]: https://github.com/MATPOWER/matpower/issues/79
[6]: https://github.com/MATPOWER/matpower/issues/90
[7]: https://osqp.org
[8]: https://github.com/MATPOWER/mptest
[9]: https://github.com/MATPOWER/mips
[10]: https://savannah.gnu.org/bugs/?52614
[11]: https://github.com/MATPOWER/most
[12]: https://github.com/ebertolazzi/mexIPOPT
[13]: https://matpower.org/doc/mpom/
[14]: https://www.artelys.com/solvers/knitro/
[15]: https://highs.dev
[16]: https://github.com/savyasachi/HiGHSMEX
[17]: https://github.com/MATPOWER/matpower/issues/15
[18]: https://github.com/MATPOWER/matpower/issues/16
