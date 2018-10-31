What's New in MIPS 1.3
----------------------

#### Released Oct 30, 2018

Below is a summary of the changes since version 1.2.2 of MIPS. See the
[`CHANGES.md`][1] file for all the gory details. For release notes for
previous versions, see Appendix C of the [MIPS User's Manual][2].

#### New Features:
  - Support for PARDISO 6.x.
  - New `mplinsolve` solver option `'LU'` for explicit LU decomposition
    with back substitution, with options in `opt.lu` for specifying the
    number of output arguments in call to `lu` (`opt.lu.nout`), whether
    to use permutation vectors or matrices (`opt.lu.vec`) and pivot
    threshold options (`opt.lu.thresh`). The following values for the
    `solver` argument act as shortcuts for specifying various
    combinations of options: `'LU3'`, `'LU3a'`, `'LU4'`, `'LU5'`,
    `'LU3m'`, `'LU3am'`, `'LU4m'`, `'LU5m'`.
    See `help mplinsolve` for details.
    *Thanks to Jose Luis Marin.*

#### Bugs Fixed:
  - Fix bug preventing `pardiso.dparm` options from being set.

#### Other Changes:
  - LaTeX source code for [MIPS User's Manual][2] included in `docs/src`.
  - Move `mplinsolve` PARDISO options to `opt.pardiso` in preparation
    for adding options for other solvers.


[1]: ../../CHANGES.md
[2]: ../MIPS-manual.pdf
