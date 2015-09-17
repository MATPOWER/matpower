========================================================
 MATPOWER - A MATLAB(R) Power System Simulation Package
========================================================

Version:    5.1

Home Page:  http://www.pserc.cornell.edu/matpower/

Authors:    Ray Zimmerman               <rz10@cornell.edu>
            Carlos E. Murillo-Sanchez   <carlos_murillo@ieee.org>
            and others, see AUTHORS file

            Fri, Mar 20, 2015

$Id$
Copyright (c) 1997-2015 by Power System Engineering Research Center (PSERC)
See http://www.pserc.cornell.edu/matpower/ for more info.

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved. This file is offered as-is,
without any warranty.


--------------
 INTRODUCTION
--------------

MATPOWER is a package of MATLAB(R) M-files for solving power flow and
optimal power flow problems. It is intended as a simulation tool for
researchers and educators that is easy to use and modify. MATPOWER
is designed to give the best performance possible while keeping the code
simple to understand and modify. It was initially developed as part
of the PowerWeb project <http://www.pserc.cornell.edu/powerweb/>.

MATPOWER can be downloaded from the MATPOWER home page above. 


--------------
 TERMS OF USE
--------------

Please see the LICENSE file for the details. But here is the summary:

- Beginning with version 5.1, the code in MATPOWER is distributed under
  the 3-clause BSD license.

- MATPOWER is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.

- While not required by the terms of the license, we do request that
  publications derived from the use of MATPOWER explicitly acknowledge
  that fact by citing [1].

[1]  R. D. Zimmerman, C. E. Murillo-Sanchez, and R. J. Thomas,
     "MATPOWER: Steady-State Operations, Planning and Analysis Tools
     for Power Systems Research and Education," Power Systems, IEEE
     Transactions on, vol. 26, no. 1, pp. 12–19, Feb. 2011.

Note:  Versions 4.0 through 5.0 were licensed under the GPL and versions
       prior to MATPOWER 4 use a different license.


-----------------
 GETTING STARTED
-----------------

System Requirements
-------------------
    - MATLAB(R) version 7 (R14) or later
      (available from The MathWorks, Inc. http://www.mathworks.com/), or
    - GNU Octave version 3.4 or later (free software available from
      http://www.gnu.org/software/octave/).

Installation
------------
1.  Follow the download instructions on the MATPOWER home page. You
    should end up with a file named matpowerXXX.zip, where XXX depends
    on the version of MATPOWER.

2.  Unzip the downloaded file. Move the resulting matpowerXXX directory
    to the location of your choice. These files should not need to be
    modified, so it is recommended that they be kept separate from your
    own code. Let $MATPOWER denote the path to this directory.

3.  Add the following directories to your MATLAB path:
      $MATPOWER   - core MATPOWER functions
      $MATPOWER/t - test scripts for MATPOWER
      (optional) subdirectories of $MATPOWER/extras -
            additional functionality and contributed code

4.  At the MATLAB prompt, type 'test_matpower' (without the quotes) to
    run the test suite and verify that MATPOWER is properly installed
    and functioning.

Running MATPOWER
----------------
To run a simple Newton power flow on the 9-bus system specified in the
file case9.m, with the default algorithm options, at the MATLAB prompt,
type:

    runpf('case9')

To load the 30-bus system data from case30.m, increase its real power
demand at bus 2 to 30 MW, then run an AC optimal power flow with
default options, type:

    define_constants;
    mpc = loadcase('case30');
    mpc.bus(2, PD) = 30;
    runopf(mpc);

By default, the results of the simulation are pretty-printed to the
screen, but the solution can also be optionally returned in a 'results'
struct. The following example shows how simple it is, after running a DC
OPF on the 118-bus system in case118.m, to access the final objective
function value, the real power output of generator 6 and the power flow
in branch 51.

    results = rundcopf('case118');
    final_objective = results.f;
    gen6_output     = results.gen(6, PG);
    branch51_flow   = results.branch(51, PF);

For additional info, see the User's Manual and the on-line help
documentation for the various MATPOWER functions. For example:
    help runpf
    help runopf
    help mpoption
    help caseformat


---------------------------
 WHAT'S NEW IN VERSION 5.1
---------------------------

Below is a summary of the changes since version 5.0 of MATPOWER. See the
CHANGES file in the docs directory for all the gory details.

* New license:
  - Switched to the more permissive 3-clause BSD license from the
    previously used GNU General Public License (GPL) v3.0.

* New case files:
  - Added four new case files, ranging from 89 up to 9421 buses,
    representing parts of the European high voltage transmission network,
    stemming from the Pan European Grid Advanced Simulation and State
    Estimation (PEGASE) project. Thanks to Cedric Josz and colleagues from
    the French Transmission System Operator.

* New documentation:
  - Added an online function reference to the web site at
    http://www.pserc.cornell.edu/matpower/docs/ref/.

* New features:
  - Added support for using PARDISO (http://www.pardiso-project.org/)
    as linear solver for computing interior-point update steps in MIPS,
    resulting in dramatic improvements in computation time and memory use
    for very large-scale problems.
  - Added support for LP/QP solver CLP (COIN_OR Linear Programming,
    http://www.coin-or.org/projects/Clp.xml). Use 'opf.dc.solver'
    option 'CLP' or qps_clp().
  - Added support for OPTI Toolbox (http://www.i2c2.aut.ac.nz/Wiki/OPTI/)
    versions of CLP, GLPK and IPOPT solvers, providing a very simple
    installation path for some free high-performance solvers on Windows
    platforms.
  - Network reduction toolbox for creating smaller approximate network
    equivalents from a larger original case file, contributed by
    Yujia Zhu and Daniel Tylavsky.
  - Added unified interface to various solvers for mixed-integer linear
    and quadratic programming (MILP/MIQP) problems.
  - Major update to have_fcn(), which now determines and caches
    version numbers and release dates for optional packages, and includes
    ability to toggle the availability of optional functionality.
  - New and updated support for 3rd party solvers:
    - High-performance IPOPT-PARDISO solver builds from the PARDISO Project
      http://www.pardiso-project.org/index.php?p=manual (at time of release
      this is MATPOWER's highest performing solver for very large scale
      AC OPF problems)
    - OPTI Toolbox versions of CLP, GLPK, IPOPT
    - CLP
    - Gurobi 6.x
    - Knitro 9.1
    - MOSEK 7.1
    - Optimization Toolbox 7.2
        - dual-simplex algorithm for linprog()
        - intlinprog() for MILP
  - New functions:
    - mplinsolve() provides unified interface for linear system solvers,
      including PARDISO and built-in backslash operator
    - miqps_matpower() provides unified interface to multiple MILP/MIQP
      solvers.
    - miqps_clex() provides a unified MILP/MIQP interface to CPLEX.
    - miqps_glpk() provides a unified MILP interface to GLPK.
    - miqps_gurobi() provides a unified MILP/MIQP interface to Gurobi.
    - miqps_mosek() provides a unified MILP/MIQP interface to MOSEK.
    - miqps_ot() provides a unified MILP interface to intlingprog().
    - mosek_symbcon() defines symbolic constants for setting
      MOSEK options.

* Other improvements:
  - Cleaned up and improved consistency of output in printpf() for
    generation and dispatchable load constraints.
  - Modified runcpf() to gracefully handle the case when the base
    and target cases are identical (as opposed to getting lost in
    an infinite loop).
  - Optional generator and dispatchable load sections in pretty-printed
    output now include off-line units.

* Bugs fixed:
  - Fixed fatal bug in case_info() for islands with no generation.
  - Fixed fatal bug in toggle_dcline() when pretty-printing results.
    Thanks to Deep Kiran for reporting.
  - Fixed sign error on multipliers on lower bound on constraints
    in qps_clp() and qps_glpk().
  - Fixed bug in handling of interface flow limits, where multipliers
    on binding interface flow limits were off by a factor of the p.u.
    MVA base.
  - Fixed minor bug with poly2pwl(), affecting units with PMAX <= 0.
  - Fixed error in qps_mosek() in printout of selected optimizer
    when using MOSEK 7.
  - Fixed bug in hasPQcap() that resulted in ignoring generator
    capability curves for units whose reactive range increases
    as real power output increases. Thanks to Irina Boiarchuk for
    reporting.
  - Fixed several incompatibilities with Matlab versions < 7.3.


---------------
 DOCUMENTATION
---------------

There are three primary sources of documentation for MATPOWER.
    - MATLAB's 'help' command
    - MATPOWER User's Manual
    - MATPOWER Online Function Reference
      (http://www.pserc.cornell.edu/matpower/docs/ref)

The User's Manual is included in the distribution (docs/manual.pdf) or
it can be downloaded separately from
http://www.pserc.cornell.edu/matpower/manual.pdf. Previous versions are
available at http://www.pserc.cornell.edu/matpower/docs/.

Each M-file has its own documentation which can be accessed by typing at
the MATLAB prompt:

    help <name of M-file>

Documentation for the case data file format can be found by typing:

    help caseformat

If something is still unclear after checking the manual and the help,
the source code *is* the documentation. ;-)

TECH NOTES

There are also two MATPOWER Technical Notes that may be of interest:

[TN1]  R. D. Zimmerman, "Uniform Price Auctions and Optimal Power Flow,
       MATPOWER Technical Note 1, February 2010.
          http://www.pserc.cornell.edu/matpower/TN1-OPF-Auctions.pdf

[TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
       their Derivatives using Complex Matrix Notation", MATPOWER
       Technical Note 2, February 2010.
          http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf

PUBLICATIONS & PRESENTATIONS

[1]  R. D. Zimmerman, C. E. Murillo-Sanchez, and R. J. Thomas,
     "MATPOWER's Extensible Optimal Power Flow Architecture," Power
     and Energy Society General Meeting, 2009 IEEE, pp. 1-7,
     July 26-30 2009.
        http://www.pserc.cornell.edu/matpower/MATPOWER-OPF.pdf
        http://dx.doi.org/10.1109/PES.2009.5275967
     slides of presentation:
       http://www.pserc.cornell.edu/matpower/MATPOWER-OPF-slides.pdf

[2]  R. D. Zimmerman, C. E. Murillo-Sanchez, and R. J. Thomas,
     "MATPOWER: Steady-State Operations, Planning and Analysis Tools
     for Power Systems Research and Education," accepted to IEEE
     Transactions on Power Systems.
        http://www.pserc.cornell.edu/matpower/MATPOWER-paper.pdf


-------------------
 OPTIONAL PACKAGES
-------------------

There are three optional packages to enhance the performance of MATPOWER
that may be downloaded separately. MINOPF and BPMPDMEX have more
restrictive licenses than MATPOWER. Please see the individual
Terms of Use for details.

 - BPMPD_MEX   MEX-file version of the high performance BPMPD interior
               point LP and QP solver. Speeds up DC and LP-based OPF
               solvers, and improves robustness of MINOPF.
               See http://www.pserc.cornell.edu/bpmpd/

 - CLP         COIN-OR Linear Programming solver implements high performance
               simplex and barrier LP and QP solvers that MATPOWER can use
               for the DC OPF. Available from:
               http://www.coin-or.org/projects/Clp.xml.

 - CPLEX       Includes high-performance, large-scale LP and QP solvers
               that MATPOWER can use for the DC OPF. Requires the
               Matlab interface to CPLEX, available from
               http://www.ibm.com/software/integration/optimization/cplex-optimizer/.

 - GLPK        GNU Linear Programming Kit includes large-scale LP solvers
               that MATPOWER can use for the DC OPF. Available from
               http://www.gnu.org/software/glpk/ and included with Octave.

 - GUROBI      Includes high-performance, large-scale LP and QP solvers
               that MATPOWER can use for the DC OPF. Requires the
               Gurobi MEX Matlab interface, available from
               http://www.convexoptimization.com/wikimization/index.php/Gurobi_mex.

 - IPOPT       An interior point optimizer for large scale non-linear
               optimization that MATPOWER can use for both AC and DC
               OPF problems. Requires the Matlab MEX interface for
               IPOPT, available from
               http://www.coin-or.org/projects/Ipopt.xml.
               Pre-built MEX binaries for Windows available with
               OPTI Toolbox from http://www.i2c2.aut.ac.nz/Wiki/OPTI/,
               and high-performance IPOPT-PARDISO pre-built MEX binaries
               for Mac and Linux from the PARDISO Project at
               http://www.pardiso-project.org/index.php?p=manual.

 - KNITRO      A general purpose optimization solver specializing in
               nonlinear problems that MATPOWER can use for AC OPFs.
               Requires the Knitro libraries, available from
               http://www.ziena.com/.

 - MINOPF      A MINOS-based AC OPF solver implemented as a Fortran MEX file.
               See http://www.pserc.cornell.edu/minopf/

 - MOSEK       Includes high-performance, large-scale LP and QP solvers
               that MATPOWER can use for the DC OPF. Requires the Matlab
               interface to MOSEK, available from http://www.mosek.com/.

 - TSPOPF      A package of three AC OPF solvers implemented as C MEX files.
               Suitable for large scale problems.
               See http://www.pserc.cornell.edu/tspopf/

These packages are distributed separately since each has it's own
license agreement and terms of use.

Note: For Windows users, one of the best ways to get access to some high
performance solvers, without dealing with the pain of building the MEX
interfaces yourself, is to install the OPTI Toolbox by Jonathan Currie,
available at: http://www.i2c2.aut.ac.nz/Wiki/OPTI/. The installation is
simple and it includes pre-built MEX files for several of the above
solvers, including CLP, GLPK and IPOPT.


--------------
 MAILING LIST
--------------

An e-mail list <MATPOWER-L@cornell.edu> has been set up to facilitate
discussion of MATPOWER. Only list subscribers are permitted to post to
the list.

Feel free to use this list to discuss anything related to MATPOWER, to
ask questions about MATPOWER, or to provide feedback to the developers
of MATPOWER, such as bug reports, patches or ideas for improvements
(though we make no guarantees about if/when they might be included).

Also, if you have any of your own MATLAB power systems code that you
would like to contribute, feel free to contact us via this list about
making it available on the MATPOWER web site.

Joining the list
----------------
To join the MATPOWER mailing list, send an e-mail to
<MATPOWER-L-request@cornell.edu> with a single line with the word "join"
in the body of the message. You must send the request from the e-mail
address where you want to receive the list's messages. And be sure it is
a plain text e-mail, that is, with no formatting, font styles or HTML
code.

Sending mail to the list
------------------------
To send an e-mail to all of the subscribers of the MATPOWER mailing
list, simply address your e-mail to <MATPOWER-L@cornell.edu>. Only
subscribers are permitted to send e-mail to the list.

Leaving the list
----------------
You can unsubscribe from the list at any time by sending an e-mail to
<MATPOWER-L-request@cornell.edu> with a single line with the word
"leave" in the body of the message
