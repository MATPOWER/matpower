=====================================================
 MATPOWER - A Matlab Power System Simulation Package
=====================================================

Version:    4.0b1

Home Page:  http://www.pserc.cornell.edu/matpower/

Authors:    Ray Zimmerman               <rz10@cornell.edu>
            Carlos E. Murillo-Sanchez   <carlos_murillo@ieee.org>
            Deqiang (David) Gan         <dgan@zju.edu.cn>

            Thu, Dec 24, 2009

$Id$
Copyright (c) 1997-2007 by Power System Engineering Research Center (PSERC)
See http://www.pserc.cornell.edu/matpower/ for more info.


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

- MATPOWER is free of charge. Anyone may use it.
- We make no warranties, express or implied. Specifically, we make no
  guarantees regarding the correctness MATPOWER's code or its fitness for any
  particular purpose.
- Any publications derived from the use of MATPOWER must acknowledge MATPOWER
  <http://www.pserc.cornell.edu/matpower/>.
- Anyone may modify MATPOWER for their own use as long as the original
  copyright notices remain in place.
- MATPOWER may not be redistributed without written permission.
- Modified versions of MATPOWER, or works derived from MATPOWER, may not be
  distributed without written permission.


-----------------
 GETTING STARTED
-----------------

System Requirements
-------------------
    - Matlab version 6.5 or later, available from The MathWorks
      http://www.mathworks.com/.

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
      $MATPOWER   Ð core MATPOWER functions
      $MATPOWER/t Ð test scripts for MATPOWER

4.  At the Matlab prompt, type 'test matpower' (without the quotes) to
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

The User's Manual is included in the distribution (manual.pdf in docs
directory) or it can be downloaded separately from 
http://www.pserc.cornell.edu/matpower/manual.pdf.

Each M-file has its own documentation which can be accessed by typing at
the Matlab prompt:

    help <name of M-file>

Documentation for the case data file format can be found by typing:

    help caseformat

If something is still unclear after checking the manual and the help,
the source code *is* the documentation ;-)


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
<lyris@cornell.edu> with the following line in the body of the
message, where "John Doe" is replaced by your real name.

    join MATPOWER-L "John Doe"

Sending mail to the list
------------------------
To send an e-mail to all of the subscribers of the MATPOWER mailing
list, simply address your e-mail to <MATPOWER-L@cornell.edu>. Only
subscribers are permitted to send e-mail to the list.

Leaving the list
----------------
You can unsubscribe from the list at any time by sending an e-mail to
<lyris@cornell.edu> with the following line in the body of the
message.

    leave MATPOWER-L
