=====================================================
 MATPOWER - A Matlab Power System Simulation Package
=====================================================

Version:    3.0b1+

Home Page:  http://www.pserc.cornell.edu/matpower/

Authors:    Ray Zimmerman               <rz10@cornell.edu>
            Carlos E. Murillo-Sanchez   <carlos_murillo@ieee.org>
            Deqiang (David) Gan         <dgan@zju.edu.cn>

            Wed, Aug 25, 2004

$Id$
Copyright (c) 1997-2004 by Power System Engineering Research Center (PSERC)
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
    - Matlab version 5 or later
    - Matlab Optimization Toolbox (required only for some OPF algorithms)

Installation
------------
    1. Unzip the downloaded file.
    2. Place files in a location on your Matlab path.
    3. Start up Matlab.

Running MATPOWER
----------------
To run a simple Newton power flow on the 9-bus system specified in the
file case9.m, with the default algorithm options, at the Matlab prompt,
type:

    >> runpf('case9')

To run an optimal power flow on the 30-bus system whose data is in
case30.m, with the default algorithm options, at the Matlab prompt,
type:

    >> runopf('case30')

To run an optimal power flow on the same system, but with the option for
MATPOWER to shut down (decommit) expensive generators, type:

    >> runuopf('case30')

For info on other parameters and options, try these:
    >> help runpf
    >> help runopf
    >> help runuopf
    >> help mpoption
    >> help caseformat

To run the test suite, place the files in the 't' subdirectory in your
Matlab path, and type:
    >> test_matpower


---------------------------
 WHAT'S NEW IN VERSION 3.0
---------------------------

Below is a summary of the changes since version 2.0 of MATPOWER. See the
CHANGES file in the docs directory for all the gory details.

New features:
- Compatibility with Matlab 7 and Optimization Toolbox 3.
- DC power flow and CD OPF solvers added.
- Gauss-Seidel power flow solver added.
- Support for MINOS-based OPF solver added (separate package,
  see http://www.pserc.cornell.edu/minopf/ for more details)
- Multiple generators at a single bus.
- Saving of solved cases as M-files or MAT-files.
- Loading of input data from M-files, MAT-files, or structs.
- Improved decommitment algorithm.
- Added a very incomplete test suite.
- Handling of price sensitive loads in OPF, modeled as negative
  generators with constant power factor constraint.

Bugs fixed:
- Phase shifters shifted the wrong direction.
- Minor fixes to IEEE CDF to MATPOWER format conversion
  (reported by D. Devaraj and Venkat)
- Flows on out-of-service lines were not being zeroed out.
  (reported by Ramazan Caglar)
- Reported total inter-tie flow values and area export values
  were incorrect.
- Several other bugs in solution printouts.


---------------
 DOCUMENTATION
---------------

There are two primary sources of documentation for MATPOWER.
    - Matlab's 'help' command
    - MATPOWER User's Manual

The User's Manual is included in the distribution (manual.pdf in docs
directory) or it can be downloaded separately from the link above.

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

There are two optional packages to enhance the performance of MATPOWER
that may be downloaded separately. Both have more restrictive licenses
than MATPOWER. Please see the individual Terms of Use for details.

 - MINOPF      A MINOS-based OPF solver which is typically much faster than
               MATPOWER's other OPF solvers.
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
<listproc@cornell.edu> with the following line in the body of the
message, where "John Doe" is replaced by your real name.

    subscribe MATPOWER-L John Doe

Sending mail to the list
------------------------
To send an e-mail to all of the subscribers of the MATPOWER mailing
list, simply address your e-mail to <MATPOWER-L@cornell.edu>. Only
subscribers are permitted to send e-mail to the list.

Leaving the list
----------------
You can unsubscribe from the list at any time by sending an e-mail to
<listproc@cornell.edu> with the following line in the body of the
message.

    unsubscribe MATPOWER-L
