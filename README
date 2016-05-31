========================================================
 MATPOWER - A MATLAB(R) Power System Simulation Package
========================================================

Version:    6.0b1

Home Page:  http://www.pserc.cornell.edu/matpower/

Authors:    Ray Zimmerman               <rz10@cornell.edu>
            Carlos E. Murillo-Sanchez   <carlos_murillo@ieee.org>
            and others, see AUTHORS file

            Wed, Jun 1, 2016

Copyright (c) 1997-2016 by Power System Engineering Research Center (PSERC)
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
    own code. Let <MATPOWER> denote the path to this directory.

3.  Add the following directories to your MATLAB path:
      <MATPOWER>        - core MATPOWER functions
      <MATPOWER>/t      - test scripts for MATPOWER
      <MATPOWER>/most   - core MOST functions
      <MATPOWER>/most/t - test scripts for MOST
      (optional) subdirectories of <MATPOWER>/extras -
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


-----------------------------
 WHAT'S NEW IN VERSION 6.0b1
-----------------------------

Below is a summary of the changes since version 5.1 of MATPOWER. See the
CHANGES file in the docs directory for all the gory details.

* New Case Files:
  - Added 9 new case files, 8 cases ranging from 1888 to 6515 buses
    representing the French system, and a 13,659-bus case representing
    parts of the of the European high voltage transmission network,
    stemming from the Pan European Grid Advanced Simulation and State
    Estimation (PEGASE) project. Thanks again to Cedric Josz and
    colleagues from the French Transmission System Operator.
  - Added case145.m, IEEE 145 bus, 50 generator dynamic test case from
    http://www.ee.washington.edu/research/pstca/dyn50/pg_tcadd50.htm.

* New Features:
  - MATPOWER Optimal Scheduling Tool (MOST) is a major new feature,
    implementing a full range of optimal power scheduling problems, from a
    simple as a deterministic, single period economic dispatch problem
    with no transmission constraints to as complex as a stochastic,
    security-constrained, combined unit-commitment and multiperiod OPF
    problem with locational contingency and load-following reserves,
    ramping costs and constraints, deferrable demands, lossy storage
    resources and uncertain renewable generation.
    See docs/MOST-manual.pdf for details.
  - General mechanism for applying modifications to an existing MATPOWER
    case. See apply_changes() and idx_ct().
  - Experimental foundation for handling of ZIP load models in power flow
    (Newton, fast-decoupled only), continuation power flow, and optimal
    power flow (MIPS, fmincon, Knitro, IPOPT solvers only). Currently,
    ZIP loads can only be specified on a system-wide basis using the
    experimental options 'exp.sys_wide_zip_loads.pw' and
    'exp.sys_wide_zip_loads.qw'.
 - Support for quadprog() under GNU Octave.
 - New contributed extras:
    - Plot electrically meaningful drawings of a MATPOWER case using
      plot_mpc() in extras/misc, contributed by Paul Cuffe.
    - Find the maximum loadability limit of a system via an optimal power
      flow and dispatchable loads, using maxloadlim() in extras/maxloadlim,
      contributed by Camille Hamon.
    - Create a quadratically-constrained quadratic programming (QCQP)
      representation of the AC power flow problem using using qcqp_opf()
      in extras/misc, contributed by Cedric Josz and colleagues.
  - New functions:
    - apply_changes() and idx_ct() provide a general mechanism for
      applying modifications to an existing MATPOWER case.
    - feval_w_path() evaluates a function located at a specified path,
      outside of the Matlab path.
    - mpopt2qpopt() provides a common interface for creating options
      struct for mi/qps_matpower() from a MATPOWER options struct.
  - New function options:
    - Option to call makeB(), makeBdc(), makePTDF(), scale_load(), and
      total_load() with full case struct (mpc) instead of individual data
      matrices (bus, branch, etc.).
    - total_load(), which now computes voltage-dependent load values,
      accepts the values 'bus' and 'area' as valid values for 'load_zone'
      argument.

* Other Improvements:
  - Changed default solver order for LP, QP, MILP, MIQP problems to move
    Gurobi before CPLEX and BPMPD after OT and GLPK.
  - Added some caching to mpoption() and made minor changes to
    nested_struct_copy() to greatly decrease the overhead added by
    mpoption() when running many small problems.
  - Added "Release History" section to Appendix of manual.
  - Many new tests.

* Bugs fixed:
  - Fixed bug in toggle_dclines() that resulted in fatal error when used
    with OPF with reactive power costs. Thanks to Irina Boiarchuk.
  - Fixed fatal bug in update_mupq() affecting cases where QMIN is greater
    than or equal to QC1MIN and QC2MIN (or QMAX is less than or equal to
    QC1MAX and QC2MAX) for all generators. Thanks Jose Miguel.
  - Copying a field containing a struct to a non-struct field with
    nested_struct_copy() now overwrites rather than causing a fatal error.
  - Fixed a bug in psse_convert_xfmr() where conversion of data for
    transformers with CZ=3 was done incorrectly. Thanks to Jose Marin
    and Yujia Zhu.
  - Fixed a fatal bug in psse_convert_xfmr() affecting transformers with
    CW and/or CZ equal to 1. Thanks to Matthias Resch.
  - Fixed a crash in have_fcn() caused by changes in OPTI Toolbox v2.15
    (or possibly v2.12)
  - Commented out isolated bus 10287 in case3375wp.m.
  - Added code to DC OPF to return success = 0 for cases where the matrix
    is singular (e.g. islanded system without slack).
  - Fixed problem in have_fcn() where SeDuMi was turning off and leaving
    off all warnings.

* Incompatible Change:
  - Removed fairmax() from the public interface by moving it inside uopf(),
    the only place it was used.


---------------
 DOCUMENTATION
---------------

There are four primary sources of documentation for MATPOWER.
    - MATPOWER User's Manual
    - MOST User's Manual
    - MATPOWER Online Function Reference
      (http://www.pserc.cornell.edu/matpower/docs/ref)
    - MATLAB's 'help' command

The MATPOWER and MOST User's Manuals are included in the distribution
(docs/MATPOWER-manual.pdf and docs/MOST-manual.pdf) or they can be downloaded
separately from http://www.pserc.cornell.edu/matpower/MATPOWER-manual.pdf and
http://www.pserc.cornell.edu/matpower/MOST-manual.pdf. Previous
versions are available at http://www.pserc.cornell.edu/matpower/docs/.

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
        http://dx.doi.org/10.1109/TPWRS.2010.2051168

[3]  H. Wang, C. E. Murillo-Sánchez, R. D. Zimmerman, R. J. Thomas,
     "On Computational Issues of Market-Based Optimal Power Flow,"
     Power Systems, IEEE Transactions on, vol. 22, no. 3,
     pp. 1185-1193, Aug. 2007.
        http://dx.doi.org/10.1109/TPWRS.2007.901301

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
               http://www.pardiso-project.org/.

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
