Change history for MIPS
=======================


Version 1.5.1 - *May 10, 2024*
------------------------------

#### 5/10/24
  - Release 1.5.1.

#### 4/16/24
  - Add `full_hist` option to enable saving in `output.hist` trajectories
    of `x`, `z`, `g`, `h`, `lam`, and `mu` (issue [#4][3]).

#### 3/12/24
  - Add Sphinx-based [Reference documentation][2].


Version 1.5 - *Dec 12, 2022*
----------------------------

#### 12/12/22
  - Release 1.5.

#### 10/28/22
  - Add to `mplinsolve()` the ability to return a struct containing the matrix
    LU factorization, and to reuse this pre-factored matrix to solve
    additional systems with different right-hand-sides by passing the struct
    in place of the A matrix to subsequent calls.
  - Add option to `mplinsolve()` to solve transposed systems by setting
    `opt.tr` to 1, including when providing the pre-factored matrix for
    the original, non-transposed system.


Version 1.4 - *Oct 8, 2020*
---------------------------

#### 10/8/20
  - Release 1.4.

#### 9/23/20
  - Add feature detection functions for LU and PARDISO, defining tags
    `'lu_vec'`, `'pardiso_legacy'`, `'pardiso_object'` and `'pardiso'`
    for `have_feature()`. Requires MP-Test 7.1 or later.
  - Remove custom LU and PARDISO detection code from `mplinsolve` and
    test functions, and switch to using `have_feature()`.

#### 8/20/20
  - In User's Manual, add missing `linsolver` option in Table 3-3.
  - In User's Manual, add Sections 4 Other Functions, 4.1 Linear System
    Solvers - mplinsolve, and move Section 3.3 Quadratic Programming
    Solver to 4.2.

#### 8/18/20
  - Silence inadvertent output from `mplinsolve()` when called without
    `solver` input argument.
  - Fix fatal errors when `mplinsolve()` is called with `LU` solver and
    dense `A` matrix.

#### 4/23/20
  - Remove additional `have_fcn()` dependencies in `mips()`,
    `t_mips_pardiso()` and `t_qps_mips()`.
  - **INCOMPATIBLE CHANGE**: Calling `mips()` with `opt.linsolver`
    set to `'PARDISO'` now results in a fatal error if PARDISO is not
    installed, rather than warning and continuing with the default
    linear solver.

#### 8/15/19
  - Fix bug in initialization of `mu` that resulted in fatal error
    in `t_mips` for some versions of MATLAB (e.g. R2013b on Mac).


Version 1.3.1 - *Jun 20, 2019*
------------------------------

#### 6/20/19
  - Release 1.3.1.
  - Add `CITATION` file.
  - Other miscellaneous documentation updates, e.g. MATPOWER website
    links updated to https://matpower.org.


Version 1.3 - *Oct 30, 2018*
----------------------------

#### 10/30/18
  - Release 1.3.

#### 9/10/18
  - Add support for PARDISO v6.x.

#### 5/7/18
  - Ignore warnings from Octave 4.4 about calling `lu` with less
    than 4 outputs and a sparse input.

#### 11/29/17
  - Fix typo preventing `pardiso.dparm` options from being set.
  - Fix incorrect use of PARDISO options in `t_mplinsolve()`.

#### 3/16/17
  - Move `mplinsolve` PARDISO options to `opt.pardiso` in preparation
    for adding options for other solvers.
  - Add `mplinsolve` solver option `'LU'` for explicit LU decomposition
    with back substitution, with options in `opt.lu` for specifying the
    number of output arguments in call to `lu` (`opt.lu.nout`),
    whether to use permutation vectors or matrices (`opt.lu.vec`) and
    pivot threshold options (`opt.lu.thresh`). The following values for
    the `solver` argument act as shortcuts for specifying various
    combinations of options: `'LU3'`, `'LU3a'`, `'LU4'`, `'LU5'`,
    `'LU3m'`, `'LU3am'`, `'LU4m'`, `'LU5m'`. See `help mplinsolve` for
    details. *Thanks to Jose Luis Marin.*

#### 12/21/16
  - Add [Travis-CI][1] integration. *Thanks to Richard Lincoln.*


Version 1.2.2 - *Dec 16, 2016* (released with MATPOWER 6.0)
-----------------------------------------------------------

#### 12/15/16
  - Moved development to GitHub: <https://github.com/MATPOWER/mips>.

#### 12/9/16
  - Renamed from MATLAB Interior Point Solver to MATPOWER Interior
    Point Solver.

#### 12/6/16
  - Remove dependence of `t_mpsolve()` on presence of `have_fcn()` to
    detect PARDISO installation.

#### 11/1/16
  - Released MATPOWER 6.0b2.


Version 1.2.1 - *Jun 1, 2016* (released with MATPOWER 6.0b1)
------------------------------------------------------------

#### 3/27/15
  - Fixed issue where default value of `feastol` option was not being
    set correctly in `mips()` when called directly (or via `qps_mips()`)
    with `feastol = 0`.


Version 1.2 - *Mar 20, 2015* (released with MATPOWER 5.1)
---------------------------------------------------------

#### 3/19/15
  - Added support for using PARDISO as linear solver for
    computing interior-point update steps in MIPS, via new `mplinsolver()`
    function and `linsolver` option.

#### 2/25/15
  - Switch to more permissive 3-clause BSD license from GPL 3.0.


Version 1.1 - *Dec 17, 2014* (released with MATPOWER 5.0)
---------------------------------------------------------

#### 12/2/14
  - Additional user-settable options: `xi`, `sigma`, `z0`,
    `alpha_min`, `rho_min`, `rho_max`, `mu_threshold` and `max_stepsize`.
  - **INCOMPATIBLE CHANGE**: The name of the option to `mips()`
    to specify the maximum number of step-size reductions when `step_control`
    is on was changed from `max_red` to `sc.red_it` for consistency with
    other MATPOWER options.


Version 1.0.2 - *Nov 5, 2013* (released with MATPOWER 5.0b1, *Jul 1, 2014*)
-------------------------------------------------------------------------

#### 11/5/13
  - Fixed a bug in MIPS where a near-singular matrix could produce
    an extremely large Newton step, resulting in incorrectly satisfying
    the relative feasibility criterion for successful termination.


Version 1.0.1 - *Apr 30, 2012*
------------------------------

#### 4/30/12
  - Fixed fatal bug in MIPS for unconstrained, scalar problems.
    *Thanks to Han Na Gwon.*

#### 12/14/11
  - Released MATPOWER 4.1.


Version 1.0 - *Feb 7, 2011* (released with MATPOWER 4.0)
--------------------------------------------------------

  - _no change_


Version 1.0b5 - *Dec 13, 2010* (released with MATPOWER 4.0b5)
-------------------------------------------------------------

  - _no change_


Version 1.0b4 - *May 21, 2010* (released with MATPOWER 4.0b4)
-------------------------------------------------------------

#### 5/10/10
  - Modified input args for Hessian evaluation function for MIPS.
    Requires `cost_mult` as 3rd argument.

#### 4/27/10
  - Check for NaN's in update step.


Version 1.0b3 - *Apr 19, 2010* (released with MATPOWER 4.0b3)
------------------------------------------------------------

#### 4/19/10
  - Changed licensing to GNU General Public license. See
    `LICENSE` and `COPYING` files for details.

#### 4/6/10
  - GNU Octave compatibility!  (tested with Octave 3.2.3)


Version 1.0b2 - *Mar 19, 2010* (released with MATPOWER 4.0b2)
------------------------------------------------------------

#### 3/10/10
  - Added optional input arg to `mipsver()` function to
    trigger return of entire version struct with `Name`, `Version`,
    `Release` and `Date` (similar to MATLAB's `ver()` function).
  - Massive help text update to more closely match MathWorks
    conventions; function names in ALL CAPS, See also ..., Examples, etc.

#### 3/8/10
  - Added `mipsver()`.

#### 3/1/10
  - Added a version number and printing of MIPS version lines
    to verbose output.

#### 2/25/10
  - Swapped `g` and `h` (`G` and `H`) in notation to match
    convention used in previously published stuff.

#### 1/25/10
  - Added `qps_mips()`, with calling syntax similar to
    `quadprog()` from the Optimization Toolbox. The main difference
    from the `quadprog()` API is that the constraints are specified
     as`l <= A*x <= u`, instead of `A*x <= b` and `Aeq*x == beq`.

#### 1/19/10
  - Renamed the pure-MATLAB interior point solver from PDIPM to
    MIPS (MATLAB Interior Point Solver).

#### 1/18/10
  - Changed order of input args to `pdipm()`, added option
    for single input struct (like fmincon), more documentation, all
    constraints are now optional, returns `exitflag = -1` for
    `numerically failed`, output includes `message` field, `lambda`
    only includes relevant fields. Added tests for pdipm as
    standalone solver.

#### 1/12/10
  - Added history field to the output with trajectories of obj,
    termination criterion, etc.

#### 1/10/10
  - Added acknowledgement of port from Hongye Wang's code.


Version 1.0b1 - *Dec 24, 2009* (released with MATPOWER 4.0b1)
-------------------------------------------------------------

#### 11/4/09
  - Removed unnecessary `return` statement at end of all M-files. If
    anything it should be an `end` statement, but even that is
    optional, so we just let functions get terminated by the
    end-of-file or another function declaration.

#### 6/5/09
  - Break out of algorithm is any element of `x` becomes NaN.

#### 3/24/09
  - Added step-controlled PDIPM variant.

#### 3/19/09
  - Fixed some bugs in default args.

#### 3/13/09
  - Added a pure MATLAB implementation of the PDIPM (primal-dual
    interior point method) solver.

----
[1]: https://travis-ci.org
[2]: https://matpower.org/doc/mips/
[3]: https://github.com/MATPOWER/mips/issues/4
