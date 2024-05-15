MP-Opt-Model
============

[MP-Opt-Model][1] is a package of MATLAB/Octave M-files for constructing
and solving mathematical programming and optimization problems. It provides
an easy-to-use, object-oriented interface for building and solving your
model. It also includes a unified interface for calling numerous LP, QP,
mixed-integer and nonlinear solvers, with the ability to switch solvers
simply by changing an input option.

It is based on code that was originally developed by Ray D. Zimmerman of
Cornell University as part of [MATPOWER][2].


System Requirements
-------------------

*   [MATLAB][3] version 7.9 (R2009b) or later, or
*   [GNU Octave][4] version 4 or later
*   [MP-Test][5]
*   [MATPOWER Interior Point Solver (MIPS)][6]


Installation
------------

**Note to [MATPOWER][2] users:** _MP-Opt-Model and its prerequisites, MIPS
and MP-Test, are included when you install [MATPOWER][2]. There is generally
no need to install MP-Opt-Model separately. You can skip directly to step 3
to verify._

Installation and use of MP-Opt-Model requires familiarity with the basic operation
of MATLAB or Octave, including setting up your MATLAB/Octave path.

1.  Clone the repository or download and extract the zip file of the MP-Opt-Model
    distribution from the [MP-Opt-Model project page][1] to the location of your
    choice. The files in the resulting `mp-opt-model` or `mp-opt-modelXXX` directory,
    where `XXX` depends on the version of MP-Opt-Model, should not need to be
    modified, so it is recommended that they be kept separate from your
    own code. We will use `<MPOM>` to denote the path to this directory.

2.  Add the following directories to your MATLAB or Octave path:
    *   `<MPOM>/lib`
    *   `<MPOM>/lib/t`

3.  At the MATLAB/Octave prompt, type `test_mp_opt_model` to run the test suite and
    verify that MP-Opt-Model is properly installed and functioning. (Note: The
    tests require functioning installations of both [MP-Test][5] and
    [MIPS][6]) The result should resemble the following:
```
  >> test_mp_opt_model
  t_have_fcn..............ok
  t_nested_struct_copy....ok
  t_nleqs_master..........ok (30 of 150 skipped)
  t_pnes_master...........ok
  t_qps_master............ok (100 of 432 skipped)
  t_miqps_master..........ok (78 of 318 skipped)
  t_nlps_master...........ok
  t_opt_model.............ok
  t_om_solve_leqs.........ok
  t_om_solve_nleqs........ok (36 of 195 skipped)
  t_om_solve_pne..........ok
  t_om_solve_qps..........ok (81 of 388 skipped)
  t_om_solve_miqps........ok (20 of 131 skipped)
  t_om_solve_nlps.........ok
  All tests successful (3461 passed, 345 skipped of 3806)
  Elapsed time 2.28 seconds.
```

Sample Usage
------------

Suppose we have the following constrained 4-dimensional quadratic
programming (QP) problem with two 2-dimensional variables, _y_ and _z_,
and two constraints, one equality and the other inequality, along with
lower bounds on all of the variables.

```
  min  1/2 [y; z]' * Q * [y; z]
  y,z
  
subject to:
  A1 * [y; z] =  b1
       A2 * y <= u2
            y >= ymin
            z <= zmax
```

And suppose the data for the problem is provided as follows.

```matlab
%% variable initial values
y0 = [1; 0];
z0 = [0; 1];

%% variable lower bounds
ymin = [0; 0];
zmax = [0; 2];

%% constraint data
A1 = [ 6 1 5 -4 ];  b1 = 4;
A2 = [ 4 9 ];       u2 = 2;

%% quadratic cost coefficients
Q = [ 8  1 -3 -4;
      1  4 -2 -1;
     -3 -2  5  4;
     -4 -1  4  12  ];
```

Below, we will show two approaches to construct and solve the problem.
The first method, based on the the Optimization Model class `opt_model`,
allows you to add variables, constraints and costs to the model
individually. Then `opt_model` automatically assembles and solves the
full model automatically.


```matlab
%%-----  METHOD 1  -----
%% build model
om = opt_model;
om.add_var('y', 2, y0, ymin);
om.add_var('z', 2, z0, [], zmax);
om.add_lin_constraint('lincon1', A1, b1, b1);
om.add_lin_constraint('lincon2', A2, [], u2, {'y'});
om.add_quad_cost('cost', Q, []);

%% solve model
[x, f, exitflag, output, lambda] = om.solve();
```

The second method requires you to construct the parameters for the full
problem manually, then call the solver function directly.

```matlab
%%-----  METHOD 2  -----
%% assemble model parameters manually
xmin = [ymin; -Inf(2,1)];
xmax = [ Inf(2,1); zmax];
x0 = [y0; z0];
A = [ A1; A2 0 0];
l = [ b1; -Inf ];
u = [ b1;  u2  ];

%% solve model
[x, f, exitflag, output, lambda] = qps_master(Q, [], A, l, u, xmin, xmax, x0);
```

The above examples are included in `<MPOM>/lib/t/qp_ex1.m` along with
some commands to print the results, yielding the output below for
each approach:


```
f = 1.875      exitflag = 1

             var bound shadow prices
     x     lambda.lower  lambda.upper
  0.5000      0.0000        0.0000
  0.0000      5.1250        0.0000
 -0.0000      0.0000        8.7500
 -0.2500      0.0000        0.0000

constraint shadow prices
lambda.mu_l  lambda.mu_u
  1.2500       0.0000
  0.0000       0.6250
```

An options struct can be passed to the `solve` method or the
`qps_master` function to select a specific solver, control the level of
progress output, or modify a solver's default parameters.

Both approaches can be applied to each of the types of problems that
MP-Opt-Model handles, namely, LP, QP, MILP, MIQP, NLP and nonlinear equations.

There are also examples in the test files in `<MPOM>/lib/t`, as well as in
the [`opf_setup()`][12] and [`opf_execute()`][13] functions in [MATPOWER][2].


Documentation
-------------

There are two primary sources of documentation for MP-Opt-Model.

The first is the [MP-Opt-Model User's Manual][7]. It can be found in
your MP-Opt-Model distribution at `<MPOM>/docs/MP-Opt-Model-manual.pdf`
and the latest version is always available at:
<https://github.com/MATPOWER/mp-opt-model/blob/master/docs/MP-Opt-Model-manual.pdf>.

And second is the built-in `help` command. As with the built-in
functions and toolbox routines in MATLAB and Octave, you can type `help`
followed by the name of a command or M-file to get help on that particular
function. Many of the M-files in MP-Opt-Model have such documentation and this
should be considered the main reference for the calling options for each
function, e.g.: `qps_master`, `miqps_master`, and `nlps_master`.


[Citing MP-Opt-Model][10]
-------------------------

We request that publications derived from the use of MP-Opt-Model
explicitly acknowledge that fact by citing the [MP-Opt-Model User's Manual][7].
The citation and DOI can be version-specific or general, as appropriate.
For version 4.2, use:

>   R. D. Zimmerman. *MP-Opt-Model User's Manual, Version 4.2*. 2024.
    [Online]. Available: https://matpower.org/docs/MP-Opt-Model-manual-4.2.pdf  
    doi: [10.5281/zenodo.11177079](https://doi.org/10.5281/zenodo.11177079)

For a version non-specific citation, use the following citation and DOI,
with *\<YEAR\>* replaced by the year of the most recent release:

>   R. D. Zimmerman. *MP-Opt-Model User's Manual*. *\<YEAR\>*.
    [Online]. Available: https://matpower.org/docs/MP-Opt-Model-manual.pdf  
    doi: [10.5281/zenodo.3818002][11]

A list of versions of the User's Manual with release dates and
version-specific DOI's can be found via the general DOI at
https://doi.org/10.5281/zenodo.3818002.


Contributing
------------

Please see our [contributing guidelines][8] for details on how to
contribute to the project or report issues.


License
-------

MP-Opt-Model is distributed under the [3-clause BSD license][9].

----
[1]: https://github.com/MATPOWER/mp-opt-model
[2]: https://matpower.org/
[3]: https://www.mathworks.com/
[4]: https://www.gnu.org/software/octave/
[5]: https://github.com/MATPOWER/mptest
[6]: https://github.com/MATPOWER/mips
[7]: docs/MP-Opt-Model-manual.pdf
[8]: CONTRIBUTING.md
[9]: LICENSE
[10]: CITATION
[11]: https://doi.org/10.5281/zenodo.3818002
[12]: https://github.com/MATPOWER/matpower/blob/master/lib/opf_setup.m
[13]: https://github.com/MATPOWER/matpower/blob/master/lib/opf_execute.m