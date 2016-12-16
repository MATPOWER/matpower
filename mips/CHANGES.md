Change history for MIPS
=======================

1.2.2 (2016-12-16) - released with MATPOWER 6.0
-----------------------------------------------

* 2016-12-15 - Moved development to GitHub <https://github.com/MATPOWER/mips>.

* 2016-12-09 - Renamed from MATLAB Interior Point Solver to MATPOWER Interior
  Point Solver.

* 2016-12-06 - Remove dependence of `t_mpsolve()` on presence of `have_fcn()` to
  detect PARDISO installation.

* 2016-11-01 - Released MATPOWER 6.0b2.


1.2.1 (2016-06-01) - released with MATPOWER 6.0b1
-------------------------------------------------

* 2015-03-27 - Fixed issue where default value of `feastol` option was not being
  set correctly in `mips()` when called directly (or via `qps_mips()`)
  with `feastol = 0`.


1.2 (2015-03-20) - released with MATPOWER 5.1
---------------------------------------------

* 2015-03-19 - Added support for using PARDISO as linear solver for
  computing interior-point update steps in MIPS, via new `mplinsolver()`
  function and `linsolver` option.

* 2015-02-25 - Switch to more permissive 3-clause BSD license from GPL 3.0.


1.1 (2014-12-17) - released with MATPOWER 5.0
---------------------------------------------

* 2014-12-02 - Additional user-settable options: `xi`, `sigma`, `z0`,
  `alpha_min`, `rho_min`, `rho_max`, `mu_threshold` and `max_stepsize`.

* 2014-12-02 - **INCOMPATIBLE CHANGE**: The name of the option to `mips()`
  to specify the maximum number of step-size reductions when `step_control`
  is on was changed from `max_red` to `sc.red_it` for consistency with
  other MATPOWER options.


1.0.2 (2013-11-05) - released with MATPOWER 5.0b1 (2014-07-01)
--------------------------------------------------------------

* 2013-11-05 - Fixed a bug in MIPS where a near-singular matrix could produce
  an extremely large Newton step, resulting in incorrectly satisfying
  the relative feasibility criterion for successful termination.


1.0.1 (2012-04-30)
------------------

* 2012-04-30 - Fixed fatal bug in MIPS for unconstrained, scalar problems.
  Thanks to Han Na Gwon.

* 2011-12-14 - Released MATPOWER 4.1.


1.0 (2011-02-07) - released with MATPOWER 4.0
---------------------------------------------

* _no change_


1.0b5 (2010-12-13) - released with MATPOWER 4.0b5
-------------------------------------------------

* _no change_


1.0b4 (2010-05-21) - released with MATPOWER 4.0b4
-------------------------------------------------

* 2010-05-10 - Modified input args for Hessian evaluation function for MIPS.
  Requires `cost_mult` as 3rd argument.

* 2010-04-27 - Check for NaN's in update step.


1.0b3 (2010-04-19) - released with MATPOWER 4.0b3
-------------------------------------------------

* 2010-04-19 - Changed licensing to GNU General Public license. See
  `LICENSE` and `COPYING` files for details.

* 2010-04-06 - GNU Octave compatibility!  (tested with Octave 3.2.3)


1.0b2 (2010-03-19) - released with MATPOWER 4.0b2
-------------------------------------------------

* 2010-03-10 - Added optional input arg to `mipsver()` function to
  trigger return of entire version struct with `Name`, `Version`,
  `Release` and `Date` (similar to MATLAB's `ver()` function).
* 2010-03-10 - Massive help text update to more closely match MathWorks
  conventions; function names in ALL CAPS, See also ..., Examples, etc.

* 2010-03-08 - Added mipsver().

* 2010-03-01 - Added a version number and printing of MIPS version lines
  to verbose output.

* 2010-02-25 - Swapped `g` and `h` (`G` and `H`) in notation to match
  convention used in previously published stuff.

* 2010-01-25 - Added `qps_mips()`, with calling syntax similar to
  `quadprog()` from the Optimization Toolbox. The main difference
  from the `quadprog()` API is that the constraints are specified
   as`l <= A*x <= u`, instead of `A*x <= b` and `Aeq*x == beq`.

* 2010-01-19 - Renamed the pure-MATLAB interior point solver from PDIPM to
  MIPS (MATLAB Interior Point Solver).

* 2010-01-18 - Changed order of input args to `pdipm()`, added option
  for single input struct (like fmincon), more documentation, all
  constraints are now optional, returns exitflag = -1 for 'numerically
  failed', output includes 'message' field, lambda only includes
  relevant fields. Added tests for pdipm as standalone solver.

* 2010-01-12 - Added history field to the output with trajectories of obj,
  termination criterion, etc.

* 2010-01-10 - Added acknowledgement of port from Hongye Wang's code.


1.0b1 (2009-12-24) - released with MATPOWER 4.0b1
-------------------------------------------------

* 2009-11-04 - Removed unnecessary `return` statement at end of all M-files. If
  anything it should be an `end` statement, but even that is
  optional, so we just let functions get terminated by the
  end-of-file or another function declaration.

* 2009-06-05 - Break out of algorithm is any element of `x` becomes NaN.

* 2009-03-24 - Added step-controlled PDIPM variant.

* 2009-03-19 - Fixed some bugs in default args.

* 2009-03-13 - Added a pure MATLAB implementation of the PDIPM (primal-dual
  interior point method) solver.
