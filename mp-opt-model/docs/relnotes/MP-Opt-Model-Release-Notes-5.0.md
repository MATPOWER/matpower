What's New in MP-Opt-Model 5.0
------------------------------

#### Released Jul 9, 2025

Below is a summary of the changes since version 4.2 of MP-Opt-Model. See
the [`CHANGES.md`][1] file for all the gory details. For release notes
for previous versions, see Appendix C of the [MP-Opt-Model User's
Manual][2].


#### New Features:
  - Support for quadratic constraints in `opt_model` and quadratically-constrained quadratic programming (QCQP) problems, including functions `qcqps_master()`, `qcqps_gurobi()`, `qcqps_knitro()`, `qcqps_nlps()`, and more.
  *Thanks to Wilson González Vanegas.*
  - Support for the open-source [HiGHS](https://highs.dev) solver for LP, QP and MILP problems, including functions `miqps_highs()`, `qps_highs()`, `have_feature_highs()`, `highsver()`, and `highs_options()`. Uses the [HiGHSMEX](https://github.com/savyasachi/HiGHSMEX) interface for HiGHS.
  - New `relax_integer` option for `opt_model.solve()`. Set to `true` to easily solve the LP or QP relaxation of a mixed integer LP or QP.
  - Support for [Artelys Knitro](https://www.artelys.com/solvers/knitro/) solver for LP and QP problems, including functions `qps_knitro()`, `knitrover()`, and `artelys_knitro_options()`.
  *Thanks to Wilson González Vanegas.*
  - Support for [Artelys Knitro](https://www.artelys.com/solvers/knitro/) 15.x, which required changes to the prior options handling.
  - Support for conversion between objects and structs to facilitate workarounds for Octave's inability to save and load classdef objects.
  - New classes:
    - `mp.opt_model` replaces legacy `opt_model` and `mp_idx_manager` classes with a new modeling API. *The legacy classes are retained for backward compatibility.*
    - `mp.set_manager` encapsulates `mp_idx_manager` functionality into an individual field object representing a specific set type, rather than in the container class.
    - `mp.set_manager_opt_model` is a subclass of `mp.set_manager` that handles common functionality, e.g. related to handling solution data, for all of the set manager subclasses used by `opt_model`.
    - `mp.sm_lin_constraint` - set manager class for linear constraints
    - `mp.sm_quad_constraint` - set manager class for quadratic constraints
    - `mp.sm_nln_constraint` - set manager class for nonlinear constraints
    - `mp.sm_nln_cost` - set manager class for general nonlinear costs
    - `mp.sm_quad_cost` - set manager class for quadratic costs
    - `mp.sm_quad_cost_legacy` - backward compatible set manager class for quadratic costs
    - `mp.sm_variable` - set manager class for variables
  - New functions:
    - `artelys_knitro_options()` sets options for Artelys Knitro.
    - `convert_constraint_multipliers()` replaces `convert_lin_constraint_multipliers()`.
    - `convert_quad_constraint()` converts bounded quadratic constraints to equality/inequality pairs.
    - `have_feature_highs()` feature detection function for HiGHS solver.
    - `highsver()` displays version of installed HiGHS.
    - `highs_options()` sets options for HiGHS.
    - `knitrover()` displays version of installed Artelys Knitro.
    - `miqps_highs()` provides standardized interface for using HiGHS to solve MILP problems.
    - `mpopt2qcqpopt()` creates/modifies `qcqps_master` options struct from MATPOWER options struct.
    - `mp.struct2object()` converts a struct back to the object from which it was created.
    - `qcqps_gurobi()` provides standardized interface for using Gurobi to solve QCQP problems.
    - `qcqps_knitro()` provides standardized interface for using Artelys Knitro to solve QCQP problems.
    - `qcqps_nlps()` provides standardized interface for using `nlps_master` to solve QCQP problems via `fmincon`, IPOPT, Artelys Knitro, or MIPS.
    - `qcqps_master()` provides a single wrapper function for calling any of MP-Opt-Model's QCQP solvers.
    - `qps_highs()` provides standardized interface for using HiGHS to solve LP/QP problems.
    - `qps_knitro()` provides standardized interface for using Artelys Knitro to solve LP/QP problems.


#### New Documentation:
Two live scripts illustrate the use of new features.
- `milp_example1.mlx` (in `mp-opt-model/lib/t`) illustrates the use of MP-Opt-Model and the new `mp.opt_model` class to build and solve an optimization (MILP) model.
- `qcqp_example1.mlx` (in `mp-opt-model/lib/t`) illustrates the new quadratic constraint features and two methods of building and solving a quadratically-constrained quadratic programming (QCQP) model.


#### Bugs Fixed:
  - Using `miqps_master()` with `'DEFAULT'` solver to solve an LP/QP problem without a MILP/MIQP solver available incorrectly threw a fatal error stating there was no solver available.


#### Other Changes:
  - Major refactor of `mp_idx_manager` to use new `mp.set_manager` class.
  - Major refactor of `opt_model` to use new `mp.set_manager_opt_model` subclasses:
    - `mp.sm_lin_constraint` - set manager class for linear constraints
    - `mp.sm_quad_constraint` - set manager class for quadratic constraints
    - `mp.sm_nln_constraint` - set manager class for nonlinear constraints
    - `mp.sm_nln_cost` - set manager class for general nonlinear costs
    - `mp.sm_quad_cost_legacy` - backward compatible set manager class for quadratic costs
    - `mp.sm_variable` - set manager class for variables
  - Deprecate the following `opt_model` methods in favor of methods of the individual `mp.set_manager` objects contained by the `opt_model` object:
    - `add_named_set()` -- use `mp.set_manager.add()`
    - `describe_idx()` -- use `mp.set_manager.describe_idx()`
    - `getN()` -- use `mp.set_manager.get_N()`
    - `init_indexed_name()` -- use `mp.set_manager.init_indexed_name()`
    - `set_type_idx_map()` -- use `mp.set_manager.set_type_idx_map()`
    - `add_lin_constraint()` -- use `mp.sm_lin_constraint.add()`
    - `add_nln_constraint()` -- use `mp.sm_nln_constraint.add()`
    - `add_nln_cost()` -- use `mp.sm_nln_cost.add()`
    - `add_quad_cost()` -- use `mp.sm_quad_cost.add()`
    - `add_var()` -- use `mp.sm_variable.add()`
    - `eval_lin_constraint()` -- use `mp.sm_lin_constraint.eval()`
    - `eval_nln_constraint()` -- use `mp.sm_nln_constraint.eval()`
    - `eval_nln_constraint_hess()` -- use `mp.sm_nln_constraint.eval_hess()`
    - `eval_nln_cost()` -- use `mp.sm_nln_cost.eval()`
    - `eval_quad_cost()` -- use `mp.sm_quad_cost.eval()`
    - `init_indexed_name()` -- use `mp.set_manager.init_indexed_name()`
    - `params_lin_constraint()` -- use `mp.sm_lin_constraint.params()`
    - `params_nln_constraint()` -- use `mp.sm_nln_constraint.params()`
    - `params_nln_cost()` -- use `mp.sm_nln_cost.params()`
    - `params_quad_cost()` -- use `mp.sm_quad_cost.params()`
    - `params_var()` -- use `mp.sm_variable.params()`
    - `set_params()` -- use `mp.set_manager.set_params()`
    - `varsets_cell2struct()` -- use `mp.sm_variable.varsets_cell2struct()`
    - `varsets_idx()` -- use `mp.sm_variable.varsets_idx()`
    - `varsets_len()` -- use `mp.sm_variable.varsets_len()`
    - `varsets_x()` -- use `mp.sm_variable.varsets_x()`
  - Update `mosek_options()` for MOSEK 11.x compatibility.
  - Update `miqps_<solver>()` functions to avoid changing MIP solution values in price computation stage. It was rounding integer variables, potentionally causing a small discrepancy between the objective value reported by the solver and the value obtained by computing directly from the returned solution _x_.
  - Deprecate the `convert_lin_constraint_multipliers()` in favor of `convert_constraint_multipliers()`. 


#### Incompatible Changes:
  - Parsed solution information was moved from the `soln` property of the `opt_model` object to the `soln` property of the individual child `mp.set_manager_opt_model` objects. Currently it is still available at the original location, but this is now deprecated.
  - The `knitro_opts` field of the `opt` input to `nlps_master()` and `nlps_knitro()` and the `solve()` method of `opt_model` has been redesigned. It is now a raw Artelys Knitro options struct, so the `opts`, `tol_x` and `tol_f` fields are no longer valid. For `tol_x` and `tol_f`, use `xtol` and `ftol`, and the contents of `opts` should be placed directly in the top level of the `knitro_opts` field.
  - Remove support for older versions of Artelys Knitro, including all references to `ktrlink` for pre-v9 versions. Currently supports Artelys Knitro version 13.1 and later.



[1]: ../../CHANGES.md
[2]: ../MP-Opt-Model-manual.pdf
[3]: https://github.com/MATPOWER/most
[4]: https://matpower.org/doc/mpom/
[5]: https://github.com/ebertolazzi/mexIPOPT
