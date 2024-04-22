What's New in MP-Opt-Model 2.0
------------------------------

#### Released July 8, 2020

Below is a summary of the changes since version 1.0 of MP-Opt-Model. See
the [`CHANGES.md`][1] file for all the gory details. For release notes
for previous versions, see Appendix C of the [MP-Opt-Model User's
Manual][2].


#### New Features
  - Add new `'fsolve'` tag to `have_fcn()` to check for availability of
   `fsolve()` function.
  - Add `nleqs_master()` function as unified interface for solving
    nonlinear equations, including implementations for `fsolve` and
    Newton's method in functions `nleqs_fsolve()` and `nleqs_newton()`,
    respectively.
  - Add support for nonlinear equations (NLEQ) to `opt_model`. For
    problems with only nonlinear equality constraints and no costs,
    the `problem_type()` method returns `'NLEQ'` and the `solve()`
    method calls `nleqs_master()` to solve the problem.
  - New functions:
      - `mpopt2nleqopt()` creates or modifies an options struct for
        `nleqs_master()` from a MATPOWER options struct.
      - `nleqs_fsolve()` provides implementation of unified nonlinear
        equation solver interface for `fsolve`.
      - `nleqs_master()` provides a single wrapper function for calling
        any of MP-Opt-Model's nonlinear equation solvers.
      - `nleqs_newton()` provides implementation of Newton's method solver
        with a unified nonlinear equation solver interface.
      - `opt_model.params_nln_constraint()` method returns parameters for
        a named (and optionally indexed) set of nonlinear constraints.
      - `opt_model.params_nln_cost()` method returns parameters for a
        named (and optionally indexed) set of general nonlinear costs.

#### Other Changes
  - Add to `eval_nln_constraint()` method the ability to compute constraints
    for a single named set.
  - Skip evaluation of gradient if `eval_nln_constraint()` is called with
    a single output argument.
  - Remove redundant MIPS tests from `test_mp_opt_model.m`.
  - Add tests for solving LP/QP, MILP/MIQP, NLP and NLEQ problems via
    `opt_model.solve()`.
  - Add Table 6.1 of valid `have_fcn()` input tags to User's Manual.


[1]: https://github.com/MATPOWER/mp-opt-model/blob/master/CHANGES.md
[2]: https://github.com/MATPOWER/mp-opt-model/blob/master/docs/MP-Opt-Model-manual.pdf
