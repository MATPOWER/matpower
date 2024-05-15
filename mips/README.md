MIPS - MATPOWER Interior Point Solver
=====================================

The [MATPOWER Interior Point Solver (MIPS)][1] is a package of MATLAB/Octave
M-files for solving non-linear programming problems (NLPs) using a primal
dual interior point method. MIPS is based on [code written in C language][2]
by Hongye Wang as a graduate student at Cornell University for optimal
power flow applications. It was later ported to the MATLAB/Octave language by
Ray D. Zimmerman for use in [MATPOWER][3].

System Requirements
-------------------

*   [MATLAB][4] version 7 (R14) or later, or
*   [GNU Octave][5] version 3.4 or later
*   [MP-Test][6] version 7.1 or later


Installation
------------

**Note to [MATPOWER][3] users:** _MIPS and its prerequiste, MP-Test, are
included when you install [MATPOWER][3]. There is generally no need to
install them separately. You can skip directly to step 3 to verify._

Installation and use of MIPS requires familiarity with the basic operation
of MATLAB or Octave, including setting up your MATLAB/Octave path.

1.  Clone the repository or download and extract the zip file of the MIPS
    distribution from the [MIPS project page][1] to the location of your
    choice. The files in the resulting `mips` or `mipsXXX` directory,
    where `XXX` depends on the version of MIPS, should not need to be
    modified, so it is recommended that they be kept separate from your
    own code. We will use `<MIPS>` to denote the path to this directory.

2.  Add the following directories to your MATLAB or Octave path:
    *   `<MIPS>/lib`
    *   `<MIPS>/lib/t`

3.  At the MATLAB/Octave prompt, type `test_mips` to run the test suite and
    verify that MIPS is properly installed and functioning. The result
    should resemble the following:
```
  >> test_mips
  t_mplinsolve......ok (8 of 180 skipped)
  t_mips............ok
  t_mips_pardiso....ok (60 of 60 skipped)
  t_qps_mips........ok
  All tests successful (304 passed, 68 skipped of 372)
  Elapsed time 0.08 seconds.
```

Documentation
-------------

There are two primary sources of documentation for MIPS. The first is
the [MIPS User's Manual][7], which gives an overview of the capabilities
and structure of MIPS and describes the formulations behind the code. It
can be found in your MIPS distribution at `<MIPS>/docs/MIPS-manual.pdf`
and the latest version is always available at:
<https://github.com/MATPOWER/mips/blob/master/docs/MIPS-manual.pdf>.

And second is the built-in `help` command. As with the built-in
functions and toolbox routines in MATLAB and Octave, you can type `help`
followed by the name of a command or M-file to get help on that particular
function. All of the M-files in MIPS have such documentation and this
should be considered the main reference for the calling options for each
function, namely: `mips`, `mipsver`, and `qps_mips`.


Publications
------------

1.  H. Wang, C. E. Murillo-Sánchez, R. D. Zimmerman, R. J. Thomas,
     ["On Computational Issues of Market-Based Optimal Power Flow,"][10]
     *Power Systems, IEEE Transactions on*, vol. 22, no. 3,
     pp. 1185-1193, Aug. 2007.  
     doi: [10.1109/TPWRS.2007.901301][10].

2.  H. Wang, *On the Computation and Application of Multi-period
    Security-constrained Optimal Power Flow for Real-time Electricity
    Market Operations*, Ph.D. thesis, Electrical and Computer
    Engineering, Cornell University, May 2007.


[Citing MIPS][11]
-----------------

We request that publications derived from the use of the MATPOWER
Interior Point Solver (MIPS) explicitly acknowledge that fact by citing
the following 2007 paper.

>   H. Wang, C. E. Murillo-Sánchez, R. D. Zimmerman, R. J. Thomas, "On
    Computational Issues of Market-Based Optimal Power Flow," *Power Systems,
    IEEE Transactions on*, vol. 22, no. 3, pp. 1185-1193, Aug. 2007.  
    doi: [10.1109/TPWRS.2007.901301][10]

The [MATPOWER Interior Point Solver (MIPS) User's Manual][7] should also be
cited explicitly in work that refers to or is derived from its content.
The citation and DOI can be version-specific or general, as appropriate.
For version 1.5.1, use:

>   R. D. Zimmerman, H. Wang. *MATPOWER Interior Point Solver (MIPS)
    User's Manual, Version 1.5.1*. 2024. [Online].
    Available: https://matpower.org/docs/MIPS-manual-1.5.1.pdf  
    doi: [10.5281/zenodo.11176870](https://doi.org/10.5281/zenodo.11176870)

For a version non-specific citation, use the following citation and DOI,
with *\<YEAR\>* replaced by the year of the most recent release:

>   R. D. Zimmerman, H. Wang. *MATPOWER Interior Point Solver (MIPS)
    User's Manual*. *\<YEAR\>*. [Online].
    Available: https://matpower.org/docs/MIPS-manual.pdf  
    doi: [10.5281/zenodo.3236506][12]

A list of versions of the User's Manual with release dates and
version-specific DOI's can be found via the general DOI at
https://doi.org/10.5281/zenodo.3236506.


Contributing
------------

Please see our [contributing guidelines][8] for details on how to
contribute to the project or report issues.

License
-------

MIPS is distributed under the [3-clause BSD license][9].

----
[1]: https://github.com/MATPOWER/mips
[2]: http://www.pserc.cornell.edu/tspopf/
[3]: https://github.com/MATPOWER/matpower
[4]: https://www.mathworks.com/
[5]: https://www.gnu.org/software/octave/
[6]: https://github.com/MATPOWER/mptest
[7]: docs/MIPS-manual.pdf
[8]: CONTRIBUTING.md
[9]: LICENSE
[10]: https://doi.org/10.1109/TPWRS.2007.901301
[11]: CITATION
[12]: https://doi.org/10.5281/zenodo.3236506
