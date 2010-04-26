=====================================================
 MATPOWER - A Matlab Power System Simulation Package
=====================================================

Version:    4.0b3

Home Page:  http://www.pserc.cornell.edu/matpower/

Authors:    Ray Zimmerman               <rz10@cornell.edu>
            Carlos E. Murillo-Sanchez   <carlos_murillo@ieee.org>
            Deqiang (David) Gan         <dgan@zju.edu.cn>

            Mon, Apr 19, 2010

$Id$
Copyright (c) 1997-2010 by Power System Engineering Research Center (PSERC)
See http://www.pserc.cornell.edu/matpower/ for more info.

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved. This file is offered as-is,
without any warranty.

--------------
 INTRODUCTION
--------------

MATPOWER is a package of Matlab M-files for solving power flow and
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
  clarify our intention to allow MATPOWER to interface with Matlab
  itself as well as any other M-files or MEX-files a user may have
  installed, regardless of their licensing terms.

- MATPOWER is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.

- While not required by the terms of the license, we do request that
  publications derived from the use of MATPOWER explicitly acknowledge
  that fact by citing [1].

[1]  R. D. Zimmerman, C. E. Murillo-Sanchez, and R. J. Thomas,
     "MATPOWER's Extensible Optimal Power Flow Architecture," Power
     and Energy Society General Meeting, 2009 IEEE, pp. 1-7,
     July 26-30 2009.

Note:  Versions prior to MATPOWER 4 use a different license.


-----------------
 GETTING STARTED
-----------------

System Requirements
-------------------
    - Matlab version 6.5 or later (available from The MathWorks
      http://www.mathworks.com/), or
    - GNU Octave version 3.2 or later (free software available from
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

3.  Add the following directories to your Matlab path:
      $MATPOWER   - core MATPOWER functions
      $MATPOWER/t - test scripts for MATPOWER
      (optional) subdirectories of $MATPOWER/extras -
            additional functionality and contributed code

4.  At the Matlab prompt, type 'test_matpower' (without the quotes) to
    run the test suite and verify that MATPOWER is properly installed
    and functioning.

Running MATPOWER
----------------
To run a simple Newton power flow on the 9-bus system specified in the
file case9.m, with the default algorithm options, at the Matlab prompt,
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
 WHAT'S NEW IN VERSION 4.0b3
-----------------------------

Below is a summary of the changes since version 4.0b2 of MATPOWER. See the
CHANGES file in the docs directory for all the gory details.

* New features:
  - Licensed under the GNU General Public License (GPL).
  - GNU Octave compatibility!
  - Included psse2matpower Perl script for converting PSS/E data files
    to MATPOWER case format.
  - Removed consecutive-bus number requirement for scale_load() and
    total_load() functions.
* Bugs fixed:
  - Smart market code would die for cases with non-consecutive bus numbers.


-----------------------------
 WHAT'S NEW IN VERSION 4.0b2
-----------------------------

Below is a summary of the changes since version 4.0b1 of MATPOWER. See the
CHANGES file in the docs directory for all the gory details.

* New features:
  - Renamed the pure-Matlab interior point solver from PDIPM to MIPS
    (Matlab Interior Point Solver), changed order of input args,
    included option for single struct input (like fmincon). Added
    tests for MIPS as a standalone solver.
  - Added new top-level wrapper function for MATPOWER's QP solver,
    called qps_matpower(), with calling syntax similar to quadprog()
    from the Optimization Toolbox, to replace mp_qp() and mp_lp().
  - Added qps_bpmpd(), qps_mips() and qps_ot(), with interface that
    matches qps_matpower() to handle implementation for BPMPD_MEX,
    MIPS and Optimization Toolbox solvers, respectively.
  - Added OPF algorithm code to output of OPF in results.raw.output.alg.
  - Added saving history of trajectory of obj, feascond, gradcond,
    compcond, costcond, etc. for MIPS solver. See results.raw.output.hist.
  - For @opf_model, added compute_cost() method, deprecated get_var_N(),
    get_lin_N() and get_nln_N() methods, replaced with single getN() method.
  - Added cost field to OPF results struct with final values of user-defined
    costs, by named set.
  - Incorporated significant updates to User's Manual (docs/manual.pdf).
  - Massive help text updates to more closely match MathWorks conventions.
    i.e. function names in ALL CAPS.
  - Added two new Tech Notes, available from MATPOWER home page.

* Bugs fixed:
  - Per unit bug with reserve costs and prices in toggle_reserves() in 4.0b1.
  - In 4.0b1 the DC OPF with user-supplied linear constraints could yield
    incorrect multipliers on variable bounds.
  - Auction code in extras/smartmarket in all previous versions contained a
    design error which has been fixed. Prices are now scaled instead of
    shifted when modified according to specified pricing rule (e.g. LAO, FRO,
    LAB, FRB, split-the-difference, etc.). Auctions with both real and reactive
    offers/bids must be type 0 or 5, type 1 = LAO is no longer allowed.

* INCOMPATIBLE CHANGES:
  - Renamed functions used to compute AC OPF cost, constraints and
    hessian, since they are used by more than fmincon:
        costfmin --> opf_costfcn
        consfmin --> opf_consfcn
        hessfmin --> opf_hessfcn


-----------------------------
 WHAT'S NEW IN VERSION 4.0b1
-----------------------------

Below is a summary of the changes since version 3.2 of MATPOWER. See the
CHANGES file in the docs directory for all the gory details.

* New features:
  - New high-performance default solvers for AC and DC OPF, based on a
    pure-Matlab implementation of the primal-dual interior point
    methods in the optional package TSPOPF. Suitable for large systems.
  - Extensive OPF enhancements:
    - Generalized, extensible OPF formulation applies to all solvers
      (AC and DC).
    - Improved method for modifying OPF formulation and output via
      a new user-defined callback function mechanism.
    - Option to co-optimize reserves based on fixed reserve
      requirements, implemented using new callback function mechanism.
    - Option to include interface flow limits (based on DC model flows),
      implemented using new callback function mechanism.
  - Option to return full power flow or OPF solution in a single
    'results' struct.
  - Ability to read and save generalized OPF user constraints, costs
    and variable limits as well as other user data in case struct.
  - New functions:
    - scale_load() conveniently modifies multiple loads.
    - makeLODF() computes line outage distribution factors.
    - total_load() retreives total load for the entire system, a
      specific zone or bus, with options to include fixed loads,
      dispatchable loads or both.
  - Use of 'areas' data matrix is now optional.
  - Many new tests in test suite.

* Bugs fixed:
  - Branch power flow limits could be violated when using the option
    OPF_FLOW_LIM = 1.

* INCOMPATIBLE CHANGES:
  - Input/output arguments to uopf() are now consistent with opf().
  - dAbr_dV() now gives now gives partial derivatives of the
    *squared* magnitudes of flows w.r.t. V, as opposed to the
    magnitudes.


---------------
 DOCUMENTATION
---------------

There are two primary sources of documentation for MATPOWER.
    - Matlab's 'help' command
    - MATPOWER User's Manual

The User's Manual is included in the distribution (docs/manual.pdf) or
it can be downloaded separately from
http://www.pserc.cornell.edu/matpower/manual.pdf.

Each M-file has its own documentation which can be accessed by typing at
the Matlab prompt:

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
     for Power Systems Research and Education," submitted to IEEE
     Transactions on Power Systems.
        http://www.pserc.cornell.edu/matpower/MATPOWER-paper.pdf


-------------------
 OPTIONAL PACKAGES
-------------------

There are three optional packages to enhance the performance of MATPOWER
that may be downloaded separately. MINOPF and BPMPDMEX have more
restrictive licenses than MATPOWER. Please see the individual
Terms of Use for details.

 - TSPOPF      A package of three AC OPF solvers implemented as C MEX files.
               Suitable for large scale problems.
               See http://www.pserc.cornell.edu/tspopf/

 - MINOPF      A MINOS-based AC OPF solver implemented as a Fortran MEX file.
               See http://www.pserc.cornell.edu/minopf/

 - BPMPD_MEX   MEX-file version of the high performance BPMPD interior
               point LP and QP solver. Speeds up DC and LP-based OPF solvers,
               and improves robustness of MINOPF.
               See http://www.pserc.cornell.edu/bpmpd/

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

Also, if you have any of your own Matlab power systems code that you
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
