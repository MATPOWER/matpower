========================================================
 MATPOWER - A MATLAB(R) Power System Simulation Package
========================================================

Version:    5.0

Home Page:  http://www.pserc.cornell.edu/matpower/

Authors:    Ray Zimmerman               <rz10@cornell.edu>
            Carlos E. Murillo-Sanchez   <carlos_murillo@ieee.org>
            Deqiang (David) Gan         <dgan@zju.edu.cn>

            Wed, Dec 17, 2014

$Id$
Copyright (c) 1997-2014 by Power System Engineering Research Center (PSERC)
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

- Beginning with version 4, the code in MATPOWER is distributed under
  the GNU General Public License (GPL) with an exception added to
  clarify our intention to allow MATPOWER to interface with MATLAB
  as well as any other MATLAB code or MEX-files a user may have
  installed, regardless of their licensing terms.

- MATPOWER is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.

- While not required by the terms of the license, we do request that
  publications derived from the use of MATPOWER explicitly acknowledge
  that fact by citing [1].

[1]  R. D. Zimmerman, C. E. Murillo-Sanchez, and R. J. Thomas,
     "MATPOWER: Steady-State Operations, Planning and Analysis Tools
     for Power Systems Research and Education," Power Systems, IEEE
     Transactions on, vol. 26, no. 1, pp. 12–19, Feb. 2011.

Note:  Versions prior to MATPOWER 4 use a different license.


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
 WHAT'S NEW IN VERSION 5.0
---------------------------

Below is a summary of the changes since version 4.1 of MATPOWER. See the
CHANGES file in the docs directory for all the gory details.

* New features:
  - Continuation power flow with tangent predictor and Newton
    method corrector, based on code contributed by Shrirang Abhyankar
    and Alex Flueck.
  - SDP_PF, a set of applications of a semidefinite programming
    relaxation of the power flow equations, contributed by
    Dan Molzahn (see 'extras/sdp_pf'):
    - Globally optimal AC OPF solver (under certain conditions).
    - Functions to check sufficient conditions for:
      - global optimality of OPF solution
      - insolvability of power flow equations
  - PSS/E RAW data conversion to MATPOWER case format (experimental)
    based on code contributed by Yujia Zhu.
  - Brand new extensible MATPOWER options architecture based on options
    struct instead of options vector.
  - Utility routines to check network connectivity and handle islands
    and isolated buses.
  - New extension implementing DC OPF branch flow soft limits.
    See 'help toggle_softlims' for details.
  - New and updated support for 3rd party solvers:
    - CPLEX 12.6
    - GLPK
    - Gurobi 5.x
    - Ipopt 3.11.x
    - Knitro 9.x.x
    - Optimization Toolbox 7.1
  - Numerous performance enhancements.
  - New functions:
    - runcpf() for continuation power flow.
    - case_info() for summarizing system information, including network
      connectivity.
    - extract_islands() to extract a network island into a separate
      MATPOWER case.
    - find_islands() to detect network islands.
    - @opt_model/describe_idx() to identify variable, constraint or
      cost row indices to aid in debugging.
    - margcost() for computing the marginal cost of generation.
    - psse2mpc() to convert PSS/E RAW date into MATPOWER case format.
    - get_losses() to compute branch series losses and reactive charging
      injections and derivatives as functions of bus voltages.
    - New experimental functions in 'extras/misc' for computing loss factors,
      checking feasibility of solutions, converting losses to negative bus
      injections and modifying an OPF problem to make it feasible.
  - Added case5.m, a 5-bus, 5-generator example case from Rui Bo.
  - New options:
    - scale_load() can scale corresponding gencost for dispatchable loads.
    - makeJac() can return full Jacobian instead of reduced version
      used in Newton power flow updates.
    - modcost() can accept a vector of shift/scale factors.
    - total_load() can return actual or nominal values for dispatchable
      loads.
    - runpf(), runopf(), etc. can send pretty-printed output to file
      without also sending it to the screen.
    - 'out.suppress_detail' option suppresses all output except system
      summary (on by default for large cases).
    - 'opf.init_from_mpc' option forces some solvers to use user-supplied
      starting point.
    - MIPS 1.1 includes many new user-settable options.
  - Reimplementated @opf_model class as sub-class of the new
    @opt_model class, which supports indexed named sets of
    variables, constraints and costs.
  - Many new tests in test suite.

* Bugs fixed:
  - Running a power flow for a case with DC lines but no gencost
    no longer causes an error.
  - Fixed a bug in runpf() where it was using the wrong initial
    voltage magnitude for generator buses marked as PQ buses. Power
    flow of solved case was not converging in zero iterations as
    expected.
  - Fixed fatal bug in MIPS for unconstrained, scalar problems.
    Thanks to Han Na Gwon.
  - Fixed a bug in int2ext() where converting a case to internal
    ordering before calling runpf() or runopf() could result in
    a fatal error due to mismatched number of columns in internal
    and external versions of data matrices. Thanks to Nasiruzzaman
    and Shiyang Li for reporting and detailing the issue.
  - DC OPF now correctly sets voltage magnitudes to 1 p.u.
    in results.
  - Fixed a bug in MIPS where a near-singular matrix could produce
    an extremely large Newton step, resulting in incorrectly satisfying
    the relative feasibility criterion for successful termination.
  - Improved the starting point created for Ipopt, Knitro and MIPS
    for variables that are only bounded on one side.
  - Fixed bug in savecase() where the function name mistakenly
    included the path when the FNAME input included a path.
  - Fixed bugs in runpf() related to enforcing generator reactive
    power limits when all generators violate limits or when
    the slack bus is converted to PQ.
  - Fixed crash when using Knitro to solve cases with all
    lines unconstrained.
  - Fixed memory issue resulting from nested om fields when
    repeatedly running an OPF using the results of a previous
    OPF as input. Thanks to Carlos Murillo-Sanchez.
  - Fixed fatal error when uopf() shuts down all gens
    attempting to satisfy Pmin limits.
  - Reactive power output of multiple generators at a PQ bus
    no longer get re-allocated when running a power flow.
  - Fixed a bug in savecase() where a gencost matrix with extra
    columns of zeros resulted in a corrupted MATPOWER case file.
  - Fixed bug in runpf() that caused a crash for cases with
    'pf.enforce_q_lims' turned on and exactly two Q limit violations,
    one Qmax and one Qmin. Thanks to Jose Luis Marin.

* INCOMPATIBLE CHANGES:
  - Optional packages TSPOPF and MINOPF must be updated to latest
    versions.
  - Renamed cdf2matp() to cdf2mpc() and updated the interface to be
    consistent with psse2mpc().
  - Removed 'ot_opts' field, replaced with 'linprog_opts' and
    'quadprog_opts' fields in the OPT argument to qps_matpower()
    and qps_ot().
  - The name of the mips() option used to specify the maximum number
    of step-size reductions with step_control on was changed from
    'max_red' to 'sc.red_it' for consistency with other MATPOWER options.
  - Removed 'max_it' option from qps_matpower() (and friends) args.
    Use algorithm specific options to control iteration limits.
  - Changed behavior of branch angle difference limits so that
    0 is interpreted as unbounded only if both ANGMIN and ANGMAX
    are zero.
  - In results struct returned by an OPF, the value of
    results.raw.output.alg is now a string, not an old-style numeric
    alg code.
  - Removed:
    - Support for Matlab 6.x.
    - Support for constr() and successive LP-based OPF solvers.
    - Support for Gurobi 4.x/gurobi_mex() interface.
    - extras/cpf, replaced by runcpf().
    - extras/psse2matpower, replaced by psse2mpc().


---------------
 DOCUMENTATION
---------------

There are two primary sources of documentation for MATPOWER.
    - MATLAB's 'help' command
    - MATPOWER User's Manual

The User's Manual is included in the distribution (docs/manual.pdf) or
it can be downloaded separately from
http://www.pserc.cornell.edu/matpower/manual.pdf.

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
               https://projects.coin-or.org/Clp.

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
               https://projects.coin-or.org/Ipopt/.

 - KNITRO      A general purpose optimization solver specializing in
               nonlinear problems that MATPOWER can use for AC OPFs.
               Requires the Knitro libraries, available from
               http://www.ziena.com/ and the Optimization Toolbox from
               The MathWorks.

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
