![MATPOWER][logo]

A Power System Simulation Package for MATLAB and Octave
-------------------------------------------------------

- **MATPOWER Website**          - https://matpower.org
- **MATPOWER GitHub Project**   - https://github.com/MATPOWER/matpower

MATPOWER is a package of M-files for solving power flow, continuation
power flow and optimal power flow problems using MATLAB or Octave. It
is intended as a simulation tool for researchers and educators that is
easy to use and modify. MATPOWER is designed to give the best
performance possible while keeping the code simple to understand and
modify.

MATPOWER releases can be downloaded from the [MATPOWER website][1],
and the latest stable and work-in-progress versions can always be
downloaded or cloned from the [MATPOWER GitHub project][2]. The
`master` branch should always contain a stable version.


System Requirements
-------------------

For all features, including those based on the new MP-Core:
*   [MATLAB][3] version 9.0 (R2016a) or later, or
*   [GNU Octave][4] version 6.2 or later

Legacy features only (from 7.1 and earlier) are also available on:
*   [MATLAB][3] version 7.9 (R2009b) or later, or
*   [GNU Octave][4] version 4 or later


Getting MATPOWER
----------------

You can either download an official *versioned release* or you can obtain
the *current development version*, which
we also attempt to keep stable enough for everyday use. The development
version includes new features and bug fixes added since the last
versioned release.

#### Versioned Releases

Download the ZIP file of the latest official versioned release from the
[MATPOWER website][1].
**Note:** This _does_ include the [MATPOWER Extras][7d].

#### Current Development Version

There are also two options for obtaining the most recent development version
of MATPOWER from the `master` branch on GitHub.
**Note:** This does _not_ include the [MATPOWER Extras][7d].

1. Clone the [MATPOWER repository from GitHub][2].
   *Use this option if you want to be able to easily update to the current
   development release, with the latest bug fixes and new features, using a
   simple `git pull` command, or if you want to help with testing or
   or development. This requires that you have a [Git client][5] (GUI
   or command-line) installed.*
    - From the command line:
        - `git clone https://github.com/MATPOWER/matpower.git`
    - Or, from the [MATPOWER GitHub repository page][2]:
        - Click the green **Clone or download** button, then **Open in Desktop**.

2. Download a ZIP file of the MATPOWER repository from GitHub.
   *Use this option if you need features or fixes introduced since
   the latest versioned release, but you do not have access to or
   are not ready to begin using Git (but don't be afraid to
   [give Git a try][6]).*
    - Go to the [MATPOWER GitHub repository page][2].
    - Click the green **Clone or download** button, then **Download ZIP**.

See [CONTRIBUTING.md][7] for information on how to get a local copy
of your own MATPOWER fork, if you are interesting in contributing
your own code or modifications.

#### MATPOWER Docker Image

MATPOWER is also available on [Docker Hub][7a] as the pre-packaged
[Docker][7b] image tagged [matpower/matpower][7c], providing a Linux
environment with Octave, MATPOWER, and the [MATPOWER Extras][7d]
pre-installed. See the [MATPOWER-Docker page][7e] for more details.

Docker images are provided for both versioned releases and
development versions.


Installation
------------

Installation and use of MATPOWER requires familiarity with the basic
operation of MATLAB or Octave. Make sure you follow the installation
instructions for the version of MATPOWER you are installing. The process
was simplified with an install script following version 6.0.

1.  **Get a copy of MATPOWER** as described above. Clone the repository
    or download and extract the ZIP file of the MATPOWER distribution
    and place the resulting directory in the location of your choice
    and call it anything you like. We will use `<MATPOWER>` as a
    placeholder to denote the path to this directory (the one
    containing `install_matpower.m`). The files in `<MATPOWER>` should
    not need to be modified, so it is recommended that they be kept
    separate from your own code.

2.  **Run the installer.**
    - Open MATLAB or Octave and change to the `<MATPOWER>` directory.
    - Run the installer and follow the directions to add the
      required directories to your MATLAB or Octave path, by typing:

            install_matpower

3.  **That's it.** There is no step 3.
    - But, if you chose not to have the installer run the test suite for
      you in step 2, you can run it now to verify that MATPOWER is
      installed and functioning properly, by typing:

            test_matpower


Running MATPOWER
----------------
To run a simple Newton power flow on the 9-bus system specified in
the file `case9.m`, with the default algorithm options, at the
MATLAB or Octave prompt, type:

```matlab
runpf('case9')
```

To load the 30-bus system data from `case30.m`, increase its real power
demand at bus 2 to 30 MW, then run an AC optimal power flow with
default options, type:

```matlab
define_constants;
mpc = loadcase('case30');
mpc.bus(2, PD) = 30;
runopf(mpc);
```

By default, the results of the simulation are pretty-printed to the
screen, but the solution can also be optionally returned in a `results`
struct. The following example shows how simple it is, after running a DC
OPF on the 118-bus system in `case118.m`, to access the final objective
function value, the real power output of generator 6 and the power flow
in branch 51.

```matlab
results = rundcopf('case118');
final_objective = results.f;
gen6_output     = results.gen(6, PG);
branch51_flow   = results.branch(51, PF);
```

For additional info, see the [MATPOWER User's Manual][8], the [on-line
function reference][9], or the built-in help documentation for the various
MATPOWER functions. For example:

    help runpf
    help runopf
    help mpoption
    help caseformat


Documentation
-------------

There are a number of sources of documentation for MATPOWER

#### User's Manuals

The User's Manuals are available as PDF files in the MATPOWER distribution
as well as online.
  - [MATPOWER User's Manual][8] -- `docs/MATPOWER-manual.pdf`
  - [MOST User's Manual][10] -- `most/docs/MOST-manual.pdf`
  - [MP-Opt-User's Manual][10a] -- `mp-opt-model/docs/MP-Opt-Model-manual.pdf`
  - [MIPS Manual][10b] -- `mips/docs/MIPS-manual.pdf`
  - [MP-Test README][10c] -- `mptest/README.md`

Current and past versions of the manuals are also available online at:
  - [https://matpower.org/doc/manuals][10d]

#### MATPOWER Documentation website

The new [MATPOWER Documentation site][10e] is intended to be the home for all
future MATPOWER documentation. It is very much a work-in-progress and
currently contains only the new:
  - [MATPOWER Developer's Manual][10f]

As new documentation is written and legacy manuals are rewritten, they
will be found here in HTML and PDF formats. The site is generated by
[Sphinx][10g] and the content is written in reStructuredText (reST) format.

#### [MATPOWER Online Function Reference][9]

#### Built-in Help

Each M-file has its own documentation which can be accessed by typing at
the MATLAB/Octave prompt:

    help <name of M-file>

Documentation for the case data file format can be found by typing:

    help caseformat

If something is still unclear after checking the manual and the help,
the source code *is* the documentation. :wink:

#### Changes

Changes to MATPOWER in each released version are summarized in the
[release notes](docs/relnotes), found in `docs/relnotes` and in
Appendix H of the [MATPOWER User's Manual][8]. A complete, detailed
change log, even for unreleased versions, is available in the
[`CHANGES.md`][11] file.


Contributing
------------

Please see our [contributing guidelines][7] for details on how to
contribute to the project or report issues.


Sponsoring the MATPOWER Project
-------------------------------

If you have found MATPOWER to be valuable, please consider supporting
the project by [becoming a sponsor](https://matpower.org/sponsor).
MATPOWER development and support require significant resources. Any
contributions from the community or other sponsors free us to focus on
that support and the development of valuable new features.


Publications and Tech Notes
---------------------------

1.  R. D. Zimmerman, C. E. Murillo-Sanchez, and R. J. Thomas,
    ["MATPOWER: Steady-State Operations, Planning and Analysis Tools
    for Power Systems Research and Education,"][12] *Power Systems, IEEE
    Transactions on*, vol. 26, no. 1, pp. 12–19, Feb. 2011.  
    doi: [10.1109/TPWRS.2010.2051168][13].

2.  R. D. Zimmerman, C. E. Murillo-Sanchez, and R. J. Thomas,
    ["MATPOWER's Extensible Optimal Power Flow Architecture,"][14]
    *Power and Energy Society General Meeting, 2009 IEEE*, pp. 1-7,
    July 26-30 2009.  
    doi: [10.1109/PES.2009.5275967][15].
     - [slides of presentation][16]

3.  H. Wang, C. E. Murillo-Sánchez, R. D. Zimmerman, R. J. Thomas,
    ["On Computational Issues of Market-Based Optimal Power Flow,"][17]
    *Power Systems, IEEE Transactions on*, vol. 22, no. 3,
    pp. 1185-1193, Aug. 2007.  
    doi: [10.1109/TPWRS.2007.901301][17].

4.  C. E. Murillo-Sanchez, R. D. Zimmerman, C. L. Anderson, and
    R. J. Thomas, ["Secure Planning and Operations of Systems with
    Stochastic Sources, Energy Storage and Active Demand,"][18]
    *Smart Grid, IEEE Transactions on*, vol. 4, no. 4, pp. 2220–2229,
    Dec. 2013.  
    doi: [10.1109/TSG.2013.2281001][18].

5.  A. J. Lamadrid, D. Munoz-Alvarez, C. E. Murillo-Sanchez,
    R. D. Zimmerman, H. D. Shin and R. J. Thomas, ["Using the MATPOWER
    Optimal Scheduling Tool to Test Power System Operation Methodologies
    Under Uncertainty,"][19] *Sustainable Energy, IEEE Transactions on*,
    vol. 10, no. 3, pp. 1280-1289, July 2019.
    doi: [10.1109/TSTE.2018.2865454][19].

6.  R. D. Zimmerman, ["Uniform Price Auctions and Optimal
    Power Flow,"][20] *MATPOWER Technical Note 1*, February 2010.  
    Available: https://matpower.org/docs/TN1-OPF-Auctions.pdf  
    doi: [10.5281/zenodo.3237850](https://doi.org/10.5281/zenodo.3237850).

7.  R. D. Zimmerman, ["AC Power Flows, Generalized OPF Costs
    and their Derivatives using Complex Matrix Notation,"][21]
    *MATPOWER Technical Note 2*, February 2010.  
    Available:
    https://matpower.org/docs/TN2-OPF-Derivatives.pdf  
    doi: [10.5281/zenodo.3237866](https://doi.org/10.5281/zenodo.3237866).

8.  B. Sereeter and R. D. Zimmerman, ["Addendum to AC Power Flows and
    their Derivatives using Complex Matrix Notation: Nodal Current
    Balance,"][22] *MATPOWER Technical Note 3*, April 2018.  
    Available: https://matpower.org/docs/TN3-More-OPF-Derivatives.pdf  
    doi: [10.5281/zenodo.3237900](https://doi.org/10.5281/zenodo.3237900).

9.  B. Sereeter and R. D. Zimmerman, ["AC Power Flows, Generalized
    OPF Costs and their Derivatives using Complex Matrix Notation
    and Cartesian Coordinate Voltages,"][23] *MATPOWER Technical
    Note 4*, April 2018.  
    Available:
    https://matpower.org/docs/TN4-OPF-Derivatives-Cartesian.pdf  
    doi: [10.5281/zenodo.3237909](https://doi.org/10.5281/zenodo.3237909).


[Citing MATPOWER][31]
---------------------

We request that publications derived from the use of MATPOWER, or the
included data files, explicitly acknowledge that fact by citing the
appropriate paper(s) and the software itself.

#### Papers

All publications derived from the use of MATPOWER, or the included data
files, should cite the 2011 MATPOWER paper:

>   R. D. Zimmerman, C. E. Murillo-Sanchez, and R. J. Thomas, "MATPOWER:
    Steady-State Operations, Planning and Analysis Tools for Power Systems
    Research and Education," *Power Systems, IEEE Transactions on*, vol. 26,
    no. 1, pp. 12-19, Feb. 2011.  
    doi: [10.1109/TPWRS.2010.2051168][13]

Publications derived from the use of the [MATPOWER Optimal Scheduling
Tool (MOST)][24] should cite the 2013 MOST paper, in addition to the
2011 MATPOWER paper above.

>   C. E. Murillo-Sanchez, R. D. Zimmerman, C. L. Anderson, and R. J. Thomas,
    "Secure Planning and Operations of Systems with Stochastic Sources,
    Energy Storage and Active Demand," *Smart Grid, IEEE Transactions on*,
    vol. 4, no. 4, pp. 2220-2229, Dec. 2013.  
    doi: [10.1109/TSG.2013.2281001][18]

Work making specific reference to the [MATPOWER Interior Point Solver
(MIPS)][32] should also cite:

>   H. Wang, C. E. Murillo-Sánchez, R. D. Zimmerman, R. J. Thomas, "On
    Computational Issues of Market-Based Optimal Power Flow," *Power Systems,
    IEEE Transactions on*, vol. 22, no. 3, pp. 1185-1193, Aug. 2007.  
    doi: [10.1109/TPWRS.2007.901301][17]

NOTE: Some of the case files included with MATPOWER request the citation
of additional publications. This includes the ACTIVSg, PEGASE, and RTE
cases. Details are available in the help text at the top of the
corresponding case files.

#### Software

For the sake of reproducibility of research results, it is best to cite
the specific version of the software used, with the version-specfic DOI.
For example, for version 8.0 of MATPOWER, use:

>   R. D. Zimmerman, C. E. Murillo-Sanchez (2024). *MATPOWER (Version 8.0)*
    [Software]. Available: https://matpower.org  
    doi: [10.5281/zenodo.11212330](https://doi.org/10.5281/zenodo.11212330)

To cite the MATPOWER software generally, without reference to a specific
version, use the following citation and DOI, with *\<YEAR\>* replaced by the
year of the most recent release:

>   R. D. Zimmerman, C. E. Murillo-Sanchez (*\<YEAR\>*). *MATPOWER*
    [Software]. Available: https://matpower.org  
    doi: [10.5281/zenodo.3236535][33]

A list of versions with release dates and version-specific DOI's can be
found via the general DOI at https://doi.org/10.5281/zenodo.3236535.

#### User's Manuals

The MATPOWER, MIPS and MOST User's Manuals should also be cited
explicitly in work that refers to or is derived from their content. As
with the software, the citation and DOI can be version-specific or
general, as appropriate. For version 8.0 of the [MATPOWER User's Manual][8],
use:

>   R. D. Zimmerman, C. E. Murillo-Sanchez. *MATPOWER User's Manual,
    Version 8.0.* 2024.  
    [Online]. Available: https://matpower.org/docs/MATPOWER-manual-8.0.pdf  
    doi: [10.5281/zenodo.11212313](https://doi.org/10.5281/zenodo.11212313)

For a version non-specific citation, use the following citation and DOI,
with *\<YEAR\>* replaced by the year of the most recent release:

>   R. D. Zimmerman, C. E. Murillo-Sanchez. *MATPOWER User's Manual.* *\<YEAR\>*.  
    [Online]. Available: https://matpower.org/docs/MATPOWER-manual.pdf  
    doi: [10.5281/zenodo.3236519][34]

A list of versions of the User's Manual with release dates and
version-specific DOI's can be found via the general DOI at
https://doi.org/10.5281/zenodo.3236519.

For information on citing the MIPS or MOST User's Manuals, please see
the [`mips/CITATION`][35] and [`most/CITATION`][36] files, respectively.

#### Recommendation

In the interest of facilitating research reproducibility and thereby
increasing the value of your MATPOWER-related research publications, we
strongly encourage you to also publish, whenever possible, all of the
code and data required to generate the results you are publishing.
[Zenodo/GitHub][37] and [IEEE DataPort][38] are two of [many available
options][39].


E-mail Lists
------------

There are two e-mail lists available to serve the MATPOWER community:

- [**Discussion List**][26] ([MATPOWER-L][26]) – to facilitate discussion
  among MATPOWER users and provide a forum for help with MATPOWER
  related questions

- [**Developer List**][27] ([MATPOWER-DEV-L][27]) – to provide a forum
  for discussion among MATPOWER users and developers related to the
  development of the MATPOWER software or proposed contributions

For details see the [Mailing Lists section][28] of the
[MATPOWER website][1].

Please select the most appropriate list for your post and do *not*
cross-post to both Discussion and Developer lists. Bug reports,
software patches, proposed enhancements, etc. should be submitted to
the [issue tracker on GitHub][29].


Optional Packages
-----------------

There are numerous optional packages to enhance the performance of
MATPOWER that must be installed separately. The terms of use and
license agreements vary. Some are free of charge for all to use,
others are only free for academic use, and others may require a
commercial license. Please see Appendix G of the [MATPOWER User's
Manual][8] for details.


License and Terms of Use
------------------------

MATPOWER is distributed as open-source under the [3-clause BSD license][30].

---

[1]: https://matpower.org
[2]: https://github.com/MATPOWER/matpower
[3]: https://www.mathworks.com/
[4]: https://www.gnu.org/software/octave/
[5]: https://git-scm.com/downloads
[6]: https://git-scm.com
[7]: CONTRIBUTING.md
[7a]: https://hub.docker.com/
[7b]: https://www.docker.com
[7c]: https://hub.docker.com/r/matpower/matpower
[7d]: https://github.com/MATPOWER/matpower-extras
[7e]: docker/MATPOWER-Docker.md
[8]: docs/MATPOWER-manual.pdf
[9]: https://matpower.org/docs/ref/
[10]: most/docs/MOST-manual.pdf
[10a]: mp-opt-model/docs/MP-Opt-Model-manual.pdf
[10b]: mips/docs/MIPS-manual.pdf
[10c]: mptest/README.md
[10d]: https://matpower.org/doc/manuals
[10e]: https://matpower.org/documentation/
[10f]: https://matpower.org/documentation/dev-manual/
[10g]: https://www.sphinx-doc.org/
[11]: CHANGES.md
[12]: https://matpower.org/docs/MATPOWER-paper.pdf
[13]: https://doi.org/10.1109/TPWRS.2010.2051168
[14]: https://matpower.org/docs/MATPOWER-OPF.pdf
[15]: https://doi.org/10.1109/PES.2009.5275967
[16]: https://matpower.org/docs/MATPOWER-OPF-slides.pdf
[17]: https://doi.org/10.1109/TPWRS.2007.901301
[18]: https://doi.org/10.1109/TSG.2013.2281001
[19]: https://doi.org/10.1109/TSTE.2018.2865454
[20]: https://matpower.org/docs/TN1-OPF-Auctions.pdf
[21]: https://matpower.org/docs/TN2-OPF-Derivatives.pdf
[22]: https://matpower.org/docs/TN3-More-OPF-Derivatives.pdf
[23]: https://matpower.org/docs/TN4-OPF-Derivatives-Cartesian.pdf
[24]: https://github.com/MATPOWER/most
[26]: https://matpower.org/mailing-lists/#discusslist
[27]: https://matpower.org/mailing-lists/#devlist
[28]: https://matpower.org/mailing-lists
[29]: https://github.com/MATPOWER/matpower/issues
[30]: LICENSE
[31]: CITATION
[32]: https://github.com/MATPOWER/mips
[33]: https://doi.org/10.5281/zenodo.3236535
[34]: https://doi.org/10.5281/zenodo.3236519
[35]: mips/CITATION
[36]: most/CITATION
[37]: https://guides.github.com/activities/citable-code/
[38]: https://ieee-dataport.org
[39]: https://www.re3data.org

[logo]: docs/src/images/MATPOWER-md.png
