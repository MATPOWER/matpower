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
of [MATPOWER][2], a package of Matlab language M-files for solving
power flow and optimal power flow problems.

System Requirements
-------------------

*   [Matlab][3] version 7 (R14) or later, or
*   [GNU Octave][4] version 3.4 or later
*   [MATPOWER][5] version 6 or later,
*   [MP-Test][6], for running the MOST test suite (included with [MATPOWER][5])
*   A good LP/MILP, QP/MIQP solver, such as Gurobi, CPLEX, MOSEK, Matlab's
    Optimization Toolbox, or GLPK (included with Octave).


Installation
------------

Installation and use of MOST requires familiarity with the basic operation
of Matlab or Octave, including setting up your Matlab path.

If you have followed the directions for installing MATPOWER found in 
the [MATPOWER User's Manual][7], then MOST should already be installed and
the appropriate paths added to your Matlab path.

To run the test suite and verify that MOST is properly installed and
functioning, at the Matlab prompt, type `test_most`. The result
should resemble the following, possibly including extra tests,
depending on the availablility of optional packages:
```matlab
>> test_most
t_most_3b_1_1_0........ok
t_most_3b_3_1_0........ok
t_most_3b_1_1_2........ok
t_most_3b_3_1_2........ok
t_most_30b_1_1_0.......ok
t_most_30b_3_1_0.......ok
t_most_30b_1_1_17......ok
t_most_30b_3_1_17......ok
t_most_fixed_res.......ok
t_most_w_ds............ok
t_most_30b_1_1_0_uc....ok
t_most_sp..............ok
t_most_spuc............ok (576 of 720 skipped)
t_most_uc..............ok (208 of 260 skipped)
t_most_suc.............ok (148 of 185 skipped)
All tests successful (762 passed, 932 skipped of 1694)
Elapsed time 93.13 seconds.
```

Documentation
-------------

There are two primary sources of documentation for MOST. The first is
the [MOST User's Manual][8], which gives an overview of the capabilities
and structure of MOST and describes the problem formulation. It
can be found in your MATPOWER distribution at `<MATPOWER>/docs/MOST-manual.pdf`
and the latest version is always available at:
<https://github.com/MATPOWER/most/blob/master/docs/MOST-manual.pdf>.

And second is the built-in `help` command. As with the built-in
functions and toolbox routines in Matlab and Octave, you can type `help`
followed by the name of a command or M-file to get help on that particular
function. All of the M-files in MOST have such documentation and this
should be considered the main reference for the calling options for each
function.

Contributing
------------

Please see our [contributing guidelines][9] for details on how to
contribute to the project or report issues.

License
-------

MOST is distributed under the [3-clause BSD license][10].

----
[1]: https://github.com/MATPOWER/most
[2]: http://www.pserc.cornell.edu/matpower/
[3]: https://github.com/MATPOWER/matpower
[4]: http://www.mathworks.com/
[5]: https://www.gnu.org/software/octave/
[6]: https://github.com/MATPOWER/mptest
[7]: https://github.com/MATPOWER/matpower/blob/master/docs/MATPOWER-manual.pdf
[8]: https://github.com/MATPOWER/most/blob/master/docs/MOST-manual.pdf
[9]: CONTRIBUTING.md
[10]: LICENSE