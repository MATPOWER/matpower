What's New in MP-Opt-Model 3.0
------------------------------

#### Released October 8, 2020

Below is a summary of the changes since version 2.1 of MP-Opt-Model. See
the [`CHANGES.md`][1] file for all the gory details. For release notes
for previous versions, see Appendix C of the [MP-Opt-Model User's
Manual][2].


#### New Features
  - Support for [OSQP][3] solver for LP and QP problems (https://osqp.org).
  - Support for modifying parameters of an existing MP-Opt-Model object.
  - Support for extracting specific named/indexed variables, costs, constraint
    values and shadow prices, etc. from a solved  MP-Opt-Model object.
  - Results of the `solve()` method saved to the `soln` field of the
    MP-Opt-Model object.
  - Allow `v0`, `vl`, and `vu` inputs to `opt_model.add_var()` method, and
    `l` and `u` inputs to `opt_model.add_lin_constraint()` to be scalars
    that get expanded automatically to the appropriate vector dimension.
  - New functions:
    - `opt_model.set_params()` method modifies parameters for a given named
      set of existing variables, costs, or constraints of an MP-Opt-Model
      object.
    - `opt_model.get_soln()` method extracts solved results for a given
      named set of variables, constraints or costs.
    - `opt_model.parse_soln()` method returns a complete set of solution
      vector and shadow price values for a solved model.
    - `opt_model.eval_lin_constraint()` method computes the constraint values
      for the full set or an individual named subset of linear constraints.
    - `qps_osqp()` provides standardized interface for using [OSQP][3] to
      solve LP/QP problems
    - `osqp_options()` initializes options for [OSQP][3] solver
    - `osqpver()` returns/displays version information for [OSQP][3]
    - ... plus 29 individual feature detection functions for
      `have_feature()`, see Table A-7 in the [MP-Opt-Model User's Manual][2]
      for details.

#### Bugs Fixed:
  - Starting point supplied to `solve()` via `opt.x0` is no longer ignored
    for nonlinear equations.
  - Calling `params_var()` method with empty `idx` no longer results in
    fatal error.
  - For `opt_model`, incorrect evaluation of constant term has been fixed for
    vector valued quadratic costs with constant term supplied as a vector.

#### Other Changes
  - Simplified logic to determine whether a quadratic cost for an
    MP-Opt-Model object is vector vs. scalar valued. If the quadratic
    coefficient is supplied as a matrix, the cost is scalar varied,
    otherwise it is vector valued.
  - Deprecated `have_fcn()` and make it a simple wrapper around the new
    modular and extensible `have_feature()`, which has now been moved to
    [MP-Test][4] (https://github.com/MATPOWER/mptest).


[1]: ../../CHANGES.md
[2]: ../MP-Opt-Model-manual.pdf
[3]: https://osqp.org
[4]: https://github.com/MATPOWER/mptest
