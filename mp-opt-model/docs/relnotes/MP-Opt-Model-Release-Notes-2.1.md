What's New in MP-Opt-Model 2.1
------------------------------

#### Released August 25, 2020

Below is a summary of the changes since version 2.0 of MP-Opt-Model. See
the [`CHANGES.md`][1] file for all the gory details. For release notes
for previous versions, see Appendix C of the [MP-Opt-Model User's
Manual][2].


#### New Features
  - Fast-decoupled Newton's and Gauss-Seidel solvers for nonlinear equations.
  - New linear equation (`'LEQ'`) problem type for models with equal
    number of variables and linear equality constraints, no costs, and no
    inequality or nonlinear equality constraints. Solved via `mplinsolve()`.
  - The `solve()` method of `opt_model` can now automatically handle mixed
    systems of equations, with both linear and nonlinear equality constraints.
  - New core nonlinear equation solver function with arbitrary, user-defined
    update function, used to implement Gauss-Seidel and Newton solvers.
  - New functions:
      - `nleqs_fd_newton()` solves a nonlinear set of equations via a
        fast-decoupled Newton's method.
      - `nleqs_gauss_seidel()` solves a nonlinear set of equations via a
        Gauss-Seidel method.
      - `nleqs_core()` implements core nonlinear equation solver with
        arbitrary update function.


#### Incompatible Changes
  - In `output` return value from `nleqs_newton()`, changed the `normF`
    field of `output.hist` to `normf`, for consistency in using lowercase
    `f` everywhere.


[1]: https://github.com/MATPOWER/mp-opt-model/blob/master/CHANGES.md
[2]: https://github.com/MATPOWER/mp-opt-model/blob/master/docs/MP-Opt-Model-manual.pdf
