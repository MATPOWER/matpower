What's New in MIPS 1.4
----------------------

#### Released Oct 8, 2020

Below is a summary of the changes since version 1.3.1 of MIPS. See the
[`CHANGES.md`][1] file for all the gory details. For release notes for
previous versions, see Appendix C of the [MIPS User's Manual][2].

#### New Features:
  - Support for `have_feature()` from [MP-Test][3] to detect availability
    and version information for optional functionality. This is a modular,
    extensible replacement for `have_fcn()` from [MATPOWER][4] and
    [MP-Opt-Model][5].
  - Feature detection functions for `lu()` and PARDISO, defining tags
    `'lu_vec'`, `'pardiso_legacy'`, `'pardiso_object'` and `'pardiso'`
    for `have\_feature()'`.
  - New functions:
    - `have_feature_lu_vec` detects support for the `lu(..., 'vector')`
      syntax.
    - `have_feature_pardiso_legacy` detects support for the legacy (v5.x)
      PARDISO interface, with individual MEX files for factor, solve, etc.
    - `have_feature_pardiso_object` detects support for the object-oriented
      (v6.x and later) PARDISO interface.
    - `have_feature_pardiso` detects availability/version of PARDISO.

#### Bugs Fixed:
  - Silence inadvertent output from `mplinsolve()` when called without
    `solver` input argument.
  - Fix fatal errors when `mplinsolve()` is called with `'LU'` solver and
    dense `A` matrix.

#### Other Changes:
  - Requires MP-Test 7.1 or later.
  - Remove `have_fcn()` dependencies in `mips()`, `t_mips_pardiso()` and
    `t_qps_mips()`.

#### Incompatible Changes:
  - Calling `mips()` with `opt.linsolver` set to `'PARDISO'` now results in
   a fatal error if PARDISO is not installed, rather than warning and
   continuing with the default linear solver.


[1]: ../../CHANGES.md
[2]: ../MIPS-manual.pdf
[3]: https://github.com/MATPOWER/mptest
[4]: https://github.com/MATPOWER/matpower
[5]: https://github.com/MATPOWER/mp-opt-model
