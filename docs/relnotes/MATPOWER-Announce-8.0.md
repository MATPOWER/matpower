What's New in MATPOWER 8.0
--------------------------

#### Released May 17, 2024

Below are some of the highlights of the changes since version 7.1 of
MATPOWER. See the [full release notes][1] and the [`CHANGES.md`][2]
file for more details. For release notes for previous versions, see
Appendix H of the [MATPOWER User's Manual][3].


#### Major Redesign:

MATPOWER 8 introduces a major redesign and rewrite of all of the MATPOWER
internals in the form of the flexible, all-new MATPOWER object-oriented core
architecture (*MP-Core*) and new two user-level frameworks to access it.
*(Previously developed under the name [MP-Element][4] in a separate repository
at [https://github.com/MATPOWER/mp-element][4].)*

- *MP-Core* -- Provides unparalleled flexibility and customization capabilities
at all levels. Facilitates new modeling (e.g. unbalanced multiphase elements,
FACTS devices, etc.), new controls (e.g. optimization of transformer taps,
PAR angles, etc.), new problem formulations, and more.
- *Flexible Framework* -- Provides new top-level functions for running power
flow (PF), continuation power flow (CPF) and optimal power flow (OPF), along
with new MATPOWER Extension API for user access to the full customization
capability of MP-Core.
- *Legacy Framework* -- Allows MP-Core modeling to be used internally by
legacy functions, facilitating use of legacy test suite.

See the new [MATPOWER Developer's Manual][5], [MATPOWER Reference Manual][5a]
and [*MATPOWER Technical Note 5*][6] for details of the new architecture. The
User's Manual has not yet been updated for the flexible framework.


#### New Three-Phase and Hybrid Proof-of-Concept Examples

Prototype examples of PF, CPF, and OPF for:
- unbalanced three-phase models
- hybrid transmission (single-phase) / distribution (three-phase) models
- multiple problem formulations

**Note:** These are proof-of-concept only, with much work remaining to
define a full set of three-phase model elements and their respective 
parameters and data formats. _But they do work!_


#### New Features:

- [MIPS][7] 1.5.1 adds to `mplinsolve()` the ability to save an LU
  factorization and reuse it to solve for additional right-hand sides.
- [MP-Opt-Model][7a] 4.2 adds support for tracing solution curves of
  general parameterized nonlinear equations (PNE), providing generalized
  continuation power flow (CPF) capabilities to MP-Core and adds support
  for Optimization Toolbox 9.5 and later (MATLAB R2023a and later).
- [MOST][8] 1.3 improves speed and memory usage on certain problems, and
  adds calculation of expected temporal locational marginal
  price (TLMP), includes transitions into first period in ramping
  reserves, and more.
- [SimulinkMATPOWER][8a], included in [MATPOWER Extras][8b], enables the use
  of MATPOWER from MATLAB's [Simulink®][8c] environment.
  *Thanks to Lukas Ortmann.*
- New options:
  - New AC power flow solver based on `fsolve()` function, selected by
    setting `'pf.alg'` option to `'FSOLVE'`.
  - New Implicit Z-bus Gauss method power flow for distribution systems,
    selected by setting `pf.alg` option to `'ZG'`.
- New functions/methods:
  - `run_mp` - Top-level function for running any task (PF, CPF, OPF) with
    the new MP-Core and flexible framework.
  - `run_pf` - Wrapper around `run_mp` for running PF.
  - `run_cpf` - Wrapper around `run_mp` for running CPF.
  - `run_opf` - Wrapper around `run_mp` for running OPF.


#### New Case Files:

- Two new European transmission system cases. *Thanks to Florin Capitanescu.*
- Two new Swedish distribution system cases. *Thanks to Gabriel Malmer.*


#### New Documentation:

- [MATPOWER Documentation site][8d] -- new documentation website with MATPOWER
  manuals, reference docs, and how-to guides
- [MATPOWER Developer's Manual][5] -- describes the architecture of the
  new MP-Core and MATPOWER flexible framework
- [MATPOWER Reference Manual][5a] -- reference for MATPOWER functions and
  classes, especially those in the new MP-Core and MATPOWER flexible framework
- [*MATPOWER Technical Note 5*][6] "MP-Element: A Unified MATPOWER
  Element Model, with Corresponding Functions and Derivatives"


#### Other Improvements:

- New MATPOWER Docker image (named [`matpower/matpower`][9]) is
  based on the official GNU Octave image ([`gnuoctave/octave`][10]) and
  is available for multiple MATPOWER and Octave versions.
- Update versions of included packages:
  - MP-Test 8.0
  - MIPS 1.5.1
  - MP-Opt-Model 4.2
  - MOST 1.3
- Numerous bug fixes.


#### Incompatible Changes:

- Remove several deprecated functions, methods and options.


[1]: https://github.com/MATPOWER/matpower/blob/master/docs/relnotes/MATPOWER-Release-Notes-8.0.md
[2]: https://github.com/MATPOWER/matpower/blob/master/CHANGES.md
[3]: https://github.com/MATPOWER/matpower/blob/master/docs/MATPOWER-manual.pdf
[4]: https://github.com/MATPOWER/mp-element
[5]: https://matpower.org/doc/dev-manual/
[5a]: https://matpower.org/doc/ref-manual/
[6]: https://matpower.org/docs/TN5-MP-Element.pdf
[7]: https://github.com/MATPOWER/mips
[7a]: https://github.com/MATPOWER/mp-opt-model
[8]: https://github.com/MATPOWER/most
[8a]: https://github.com/MATPOWER/mx-simulink_matpower
[8b]: https://github.com/MATPOWER/matpower-extras
[8c]: https://www.mathworks.com/products/simulink.html
[8d]: https://matpower.org/doc/
[9]: https://hub.docker.com/r/matpower/matpower
[10]: https://hub.docker.com/r/gnuoctave/octave
