What's New in MIPS 1.5
----------------------

#### Released Dec 12, 2022

Below is a summary of the changes since version 1.4 of MIPS. See the
[`CHANGES.md`][1] file for all the gory details. For release notes for
previous versions, see Appendix C of the [MIPS User's Manual][2].

#### New Features:
  - Add to `mplinsolve()` the ability to return a struct containing the
    matrix LU factorization, and to reuse this pre-factored matrix to solve
    additional systems with different right-hand-sides by passing the
    struct in place of the A matrix to subsequent calls.
  - Add option to `mplinsolve()` to solve transposed systems by setting
    `opt.tr` to 1, including when providing the pre-factored matrix for
    the original, non-transposed system.


[1]: ../../CHANGES.md
[2]: ../MIPS-manual.pdf
