What's New in MP-Opt-Model 4.2
------------------------------

#### Released May 10, 2024

Below is a summary of the changes since version 4.1 of MP-Opt-Model. See
the [`CHANGES.md`][1] file for all the gory details. For release notes
for previous versions, see Appendix C of the [MP-Opt-Model User's
Manual][2].


#### New Features:
  - Option for `opt_model.add_lin_constraint()` to provide/store the
    transpose of the _A_ matrix instead of the original. This can
    potentially save significant memory for sparse matrices with many
    more columns than rows. E.g. storage constraints in [MOST][3]
    for 8760 hour planning horizon.
  - Add support to `nlps_master()` for calling `nlps_<my_solver>()` by
    setting `opt.alg` to `'<MY_SOLVER>'` to allow for handling custom
    NLP solvers.
  - Add support to `miqps_master()` for calling `miqps_<my_solver>()`
    by setting `opt.alg` to `'<MY_SOLVER>'` to allow for handling custom
    MILP/MIQP solvers.
  - Add to the `parse_soln()` method of `opt_model` an optional `stash`
    input argument that, if present and true, causes the parsed
    solution to be stored back in the object, as the `solve()` method
    was already doing when `opt.parse_soln` is true.
  - New Sphinx-based [Reference documentation][4].
  - New functions:
    - `convert_lin_constraint()` converts linear constraints from a
      single set of doubly-bounded inequality constraints to separate
      sets of equality and upper-bounded inequality constraints.
    - `convert_lin_constraint_multipliers()` converts multipliers on
      linear constraints from separate sets for equality and
      upper-bounded inequality constraints to those for doubly-bounded
      inequality constraints.
  - New `opt_model` methods:
      - `is_solved()` indicates whether the model has been solved.
      - `has_parsed_soln()` indicates whether a parsed solution is
        available in the model.
      - `display_soln()` display the results of a solved model,
        including values, bounds and shadow prices for variables and
        linear constraints, values and shadow prices for nonlinear
        constraints, and individual cost components.

#### Bugs Fixed:
  - Clear cached parameters after updating linear constraints or
    quadratic costs via `opt_model.set_params()`.
  - In `miqps_mosek()` the lower and upper bounds of binary variables
    got overwritten with 0 and 1, respectively, effectively relaxing
    any potentially tighter bounds provided as input.
  - Fix false positive in `have_feature_fsolve` in case where the file
    is present, but without a valid license.

#### Other Changes:
  - Update for compatibility with MATLAB R2023a (Optimization Toolbox
    9.5) and later, which removed `x0` as a valid input to `linprog`.
  - Update `have_feature_ipopt()` to recognize IPOPT MEX installations
    from Enrico Bertolazzi's [mexIPOPT][5], which include MEX files
    that have been renamed to architecture-specific names along with an
    `ipopt.m` wrapper function to call the appropriate one.
    *Thanks to Carlos Murillo-Sánchez.*

    _**Note:** While MP-Opt-Model no longer requires this, my
    recommendation is still to simply rename the MEX file to
    `ipopt.<mexext>`, with the appropriate architecture-specific
    extension, and delete the unnecessary `ipopt.m` entirely._
  - Always skip price computation stage in `miqps_<solver>()` functions
    for pure (as opposed to mixed) integer problems.
  - Add caching of aggregate output parameters in
    `opt_model.params_var()`.


[1]: ../../CHANGES.md
[2]: ../MP-Opt-Model-manual.pdf
[3]: https://github.com/MATPOWER/most
[4]: https://matpower.org/doc/mpom/
[5]: https://github.com/ebertolazzi/mexIPOPT
