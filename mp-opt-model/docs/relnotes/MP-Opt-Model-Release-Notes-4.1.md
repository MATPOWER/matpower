What's New in MP-Opt-Model 4.1
------------------------------

#### Released December 13, 2022

Below is a summary of the changes since version 4.0 of MP-Opt-Model. See
the [`CHANGES.md`][1] file for all the gory details. For release notes
for previous versions, see Appendix C of the [MP-Opt-Model User's
Manual][2].


#### New Features
  - Add support to `qps_master()` for calling `qps_<my_solver>()` by
    setting `opt.alg` to `'<MY_SOLVER>'` to allow for custom solvers.

#### Other Changes
  - Update for compatibility with Artelys Knitro 13.1 and later.
  - Add elapsed time in seconds to results of the `solve()` method of
    `opt_model`, returned in `om.soln.output.et`.
  - Add `runtime` field to `output` argument of `qps_glpk()` and
    `qps_mosek()`.


[1]: ../../CHANGES.md
[2]: ../MP-Opt-Model-manual.pdf
