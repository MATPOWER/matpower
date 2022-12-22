What's New in MATPOWER 8.0b1
----------------------------

#### Released Dec 22, 2022

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
at all levels. Facilitates new modeling (e.g. unbalanced multiphase elements, FACTS devices, etc.), new controls (e.g. optimization of transformer taps,
PAR angles, etc.), new problem formulations, and more.
- *Flexible Framework* -- Provides new top-level functions for running power
flow (PF), continuation power flow (CPF) and optimal power flow (OPF), along
with new MATPOWER Extension API for user access to the full customization
capability of MP-Core.
- *Legacy Framework* -- Allows MP-Core modeling to be used internally by
legacy functions, facilitating use of legacy test suite.

See the new [MATPOWER Developer's Manual][5] and [*MATPOWER Technical
Note 5*][6] for details of the new architecture. The User's Manual has
not yet been updated for the flexible framework.


#### New Features:

- [MIPS][7] 1.5 adds to `mplinsolve()` the ability to save an LU
  factorization and reuse it to solve for additional right-hand sides.
- [MOST][8] 1.2 adds calculation of expected temporal locational marginal
  price (TLMP), includes transitions into first period in ramping
  reserves, and more.
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

- Two new European case files. *Thanks to Florin Capitanescu.*


#### New Documentation:

- [MATPOWER Developer's Manual][5] -- describes the architecture of the
  new MP-Core and MATPOWER flexible framework
- [*MATPOWER Technical Note 5*][6] "MP-Element: A Unified MATPOWER
  Element Model, with Corresponding Functions and Derivatives"


#### Other Improvements:

- New MATPOWER Docker image (named [`matpower/matpower`][9]) is
  based on the official GNU Octave image ([`gnuoctave/octave`][10]) and
  is available for multiple MATPOWER and Octave versions.
- Update versions of included packages:
  - MP-Test 8.0b1.
  - MIPS 1.5.
  - MP-Opt-Model 4.1
  - MOST 1.2.
- Numerous bug fixes.


#### Incompatible Changes:

- Remove several deprecated functions, methods and options.


[1]: https://github.com/MATPOWER/matpower/blob/master/docs/relnotes/MATPOWER-Release-Notes-8.0.md
[2]: https://github.com/MATPOWER/matpower/blob/master/CHANGES.md
[3]: https://github.com/MATPOWER/matpower/blob/master/docs/MATPOWER-manual.pdf
[4]: https://github.com/MATPOWER/mp-element
[5]: https://matpower.org/documentation/dev-manual/
[6]: https://matpower.org/docs/TN5-MP-Element.pdf
[7]: https://github.com/MATPOWER/mips
[8]: https://github.com/MATPOWER/most
[9]: https://hub.docker.com/r/matpower/matpower
[10]: https://hub.docker.com/r/gnuoctave/octave
