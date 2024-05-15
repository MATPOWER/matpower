MOST - MATPOWER Optimal Scheduling Tool
=======================================

The [MATPOWER Optimal Scheduling Tool (MOST)][1] is framework for
solving generalized steady-state electric power scheduling problems.

MOST can be used to solve problems as simple as a deterministic,
single period economic dispatch problem with no transmission
constraints or as complex as a stochastic, security-constrained,
combined unit-commitment and multiperiod optimal power flow
problem with locational contingency and load-following reserves,
ramping costs and constraints, deferrable demands, lossy storage
resources and uncertain renewable generation.

While the problem formulation is general and incorporates a full
nonlinear AC network model, the current implementation is limited to
DC power flow modeling of the network. Some work has been done on an
AC implementation, but it is not yet ready for release.

The primary developers of MOST are Carlos E. Murillo-Sanchez and
Ray D. Zimmerman, with significant contributions from
Daniel Munoz-Alvarez and Alberto J. Lamadrid. It is built on top
of [MATPOWER][2], a package of MATLAB/Octave M-files for solving
power flow and optimal power flow problems.

System Requirements
-------------------

This version of MOST requires:
*   [MATPOWER][2] version 8.x or later, _(see MATPOWER system requirements
    for details of required versions of [MATLAB][4] or [GNU Octave][5])_
*   _(highly recommended)_ A high-performance LP/MILP, QP/MIQP solver,
   such as Gurobi, CPLEX, MOSEK, MATLAB's Optimization Toolbox, or GLPK
   _(included with Octave)_.


Installation
------------

The preferred method of installation is simply to install [MATPOWER][3],
which is a prerequisite for MOST and also includes its own copy of MOST.

If you have followed the directions for installing MATPOWER found in
the [MATPOWER User's Manual][6], then MOST should already be installed and
the appropriate paths added to your MATLAB/Octave path.

To run the test suite and verify that MOST is properly installed and
functioning, at the MATLAB/Octave prompt, type `test_most`. The result
should resemble the following, possibly including extra tests,
depending on the availablility of optional packages:
```
>> test_most
t_most_3b_1_1_0........ok
t_most_3b_3_1_0........ok
t_most_30b_1_1_0.......ok
t_most_30b_3_1_0.......ok
t_most_fixed_res.......ok
t_most_30b_1_1_0_uc....ok
t_most_mpopf...........ok
t_most_uc..............ok (272 of 340 skipped)
t_most_suc.............ok (148 of 185 skipped)
t_most_tlmp............ok
t_most_w_ds............ok
All tests successful (359 passed, 420 skipped of 779)
Elapsed time 39.45 seconds.
```

If, for some reason, you prefer to install your own copy of MOST directly
from the [MOST GitHub repository][1], simply clone the repository to the
location of your choice, where we use `<MOST>` to denote the path the
resulting `most` directory. Then add the following directories to your
MATLAB or Octave path:
 *  `<MOST>/lib`
 *  `<MOST>/lib/t`

It is important that they appear before MATPOWER in your path if you want
to use this version of MOST, rather than the one included with MATPOWER.


Documentation
-------------

There are two primary sources of documentation for MOST. The first is
the [MOST User's Manual][7], which gives an overview of the capabilities
and structure of MOST and describes the problem formulation. It
can be found in your MATPOWER distribution at `<MATPOWER>/most/docs/MOST-manual.pdf`
and the latest version is always available at:
<https://github.com/MATPOWER/most/blob/master/docs/MOST-manual.pdf>.

And second is the built-in `help` command. As with the built-in
functions and toolbox routines in MATLAB and Octave, you can type `help`
followed by the name of a command or M-file to get help on that particular
function. All of the M-files in MOST have such documentation and this
should be considered the main reference for the calling options for each
function.


Publications
------------

1.  R. D. Zimmerman, C. E. Murillo-Sanchez, and R. J. Thomas,
    ["MATPOWER: Steady-State Operations, Planning and Analysis Tools
    for Power Systems Research and Education,"][9] *Power Systems, IEEE
    Transactions on*, vol. 26, no. 1, pp. 12–19, Feb. 2011.  
    doi: [10.1109/TPWRS.2010.2051168][9].

2.  C. E. Murillo-Sanchez, R. D. Zimmerman, C. L. Anderson, and
    R. J. Thomas, ["Secure Planning and Operations of Systems with
    Stochastic Sources, Energy Storage and Active Demand,"][10]
    *Smart Grid, IEEE Transactions on*, vol. 4, no. 4, pp. 2220–2229,
    Dec. 2013.  
    doi: [10.1109/TSG.2013.2281001][10].

3.  A. J. Lamadrid, D. Munoz-Alvarez, C. E. Murillo-Sanchez,
    R. D. Zimmerman, H. D. Shin and R. J. Thomas, ["Using the MATPOWER
    Optimal Scheduling Tool to Test Power System Operation Methodologies
    Under Uncertainty,"][11] *Sustainable Energy, IEEE Transactions on*,
    vol. 10, no. 3, pp. 1280-1289, July 2019.  
    doi: [10.1109/TSTE.2018.2865454][11].


[Citing MOST][13]
-----------------

We request that publications derived from the use of the [MATPOWER
Optimal Scheduling Tool (MOST)][1] explicitly acknowledge that fact by
citing both the 2011 MATPOWER paper and the 2013 MOST paper.

>   R. D. Zimmerman, C. E. Murillo-Sanchez, and R. J. Thomas, "MATPOWER:
    Steady-State Operations, Planning and Analysis Tools for Power Systems
    Research and Education," *Power Systems, IEEE Transactions on*, vol. 26,
    no. 1, pp. 12-19, Feb. 2011.  
    doi: [10.1109/TPWRS.2010.2051168][9]

>   C. E. Murillo-Sanchez, R. D. Zimmerman, C. L. Anderson, and R. J. Thomas,
    "Secure Planning and Operations of Systems with Stochastic Sources,
    Energy Storage and Active Demand," *Smart Grid, IEEE Transactions on*,
    vol. 4, no. 4, pp. 2220-2229, Dec. 2013.  
    doi: [10.1109/TSG.2013.2281001][10]

The [MATPOWER Optimal Scheduling Tool (MOST) User's Manual][7] should
also be cited explicitly in work that refers to or is derived from its
content. The citation and DOI can be version-specific or general, as
appropriate. For version 1.3, use:

>   R. D. Zimmerman, C. E. Murillo-Sanchez. *MATPOWER Optimal Scheduling
    Tool (MOST) User's Manual, Version 1.3*. 2024. [Online].
    Available: https://matpower.org/docs/MOST-manual-1.3.pdf  
    doi: [10.5281/zenodo.11177189](https://doi.org/10.5281/zenodo.11177189)

For a version non-specific citation, use the following citation and DOI,
with *\<YEAR\>* replaced by the year of the most recent release:

>   R. D. Zimmerman, C. E. Murillo-Sanchez. *MATPOWER Optimal Scheduling
    Tool (MOST) User's Manual*. *\<YEAR\>*. [Online].
    Available: https://matpower.org/docs/MOST-manual.pdf  
    doi: [10.5281/zenodo.3236531][14]

A list of versions of the User's Manual with release dates and
version-specific DOI's can be found via the general DOI at
https://doi.org/10.5281/zenodo.3236531.


Contributing
------------

Please see our [contributing guidelines][8] for details on how to
contribute to the project or report issues.


License
-------

MOST is distributed under the [3-clause BSD license][12].

----
[1]: https://github.com/MATPOWER/most
[2]: https://matpower.org
[3]: https://github.com/MATPOWER/matpower
[4]: https://www.mathworks.com/
[5]: https://www.gnu.org/software/octave/
[6]: https://github.com/MATPOWER/matpower/blob/master/docs/MATPOWER-manual.pdf
[7]: docs/MOST-manual.pdf
[8]: CONTRIBUTING.md
[9]: https://doi.org/10.1109/TPWRS.2010.2051168
[10]: https://doi.org/10.1109/TSG.2013.2281001
[11]: https://doi.org/10.1109/TSTE.2018.2865454
[12]: LICENSE
[13]: CITATION
[14]: https://doi.org/10.5281/zenodo.3236531
