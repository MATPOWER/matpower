What's New in MATPOWER 8.0b1
----------------------------

#### Released Dec 22, 2022

Below is a summary of the changes since version 7.1 of MATPOWER. See the
[`CHANGES.md`][1] file for all the gory details. For release notes for
previous versions, see Appendix H of the [MATPOWER User's Manual][2].


#### Major Redesign:

MATPOWER 8 introduces a major redesign and rewrite of all of the MATPOWER
internals in the form of the flexible, all-new MATPOWER object-oriented core
architecture (*MP-Core*) and new two user-level frameworks to access it.
*(Previously developed under the name [MP-Element][3] in a separate repository
at [https://github.com/MATPOWER/mp-element][3].)*

- *MP-Core* -- Provides unparalleled flexibility and customization capabilities
at all levels. Built around an explicit three-layer modeling structure, with
*data*, *network* and *mathematical* model objects and a *task* object to
manage them. Facilitates new modeling (e.g. unbalanced multiphase elements, FACTS devices, etc.), new controls (e.g. optimization of transformer taps,
PAR angles, etc.), new problem formulations, and more.
- *Flexible Framework* -- Provides new top-level functions (note underscores)
`run_pf()`, `run_cpf()`, `run_opf()` for running power flow (PF), continuation
power flow (CPF) and optimal power flow (OPF), along with new MATPOWER
Extension API for user access to the full customization capability of MP-Core.
- *Legacy Framework* -- Allows MP-Core modeling to be used internally by
legacy `runpf()`, `runcpf()`, `runopf()`, etc. Facilitates use of legacy
test suite, but is restricted to legacy customization mechanisms.

See the new [MATPOWER Developer's Manual][4] and [*MATPOWER Technical
Note 5*][5] for details of the new architecture. The User's Manual has
not yet been updated for the flexible framework.

The features based on MP-Core are available under MATLAB 9.1 or Octave 6.2
or newer, where the legacy framework uses MP-Core's new modeling by default
for:
- `rundcpf` -- DC power flow
- `rundcopf` -- DC optimal power flow
- `runpf` -- AC power flow, for all except radial and hybrid Newton-Raphson
   formulations/solvers
- `runcpf` -- AC continuation power flow
- `runopf` -- AC OPF, for solvers MIPS, `fmincon`, IPOPT, and Artelys
   Knitro, for all formulations

Under older versions of MATLAB or Octave, MATPOWER automatically reverts to the legacy core code, which is still included. The legacy core can also be selected manually on newer versions with the `'exp.use_legacy_core'` option or by disabling MP-Core in the legacy framework via `have_feature('mp_core', 0)`.


#### New Features:

- [MP-Test][6] 8.0b1 includes new functions for testing string values
  and text file contents. For details, see the [MP-Test 8.0b1 release
  notes][7].
- [MIPS][8] 1.5 adds to `mplinsolve()` the ability to save an LU
  factorization and reuse it to solve for additional right-hand sides.
  For details, see the [MIPS 1.5 release notes][9]
  *(also in Appendix C in the [MIPS User's Manual][10])*.
- [MP-Opt-Model][11] 4.1 adds support for Artelys Knitro 13.1 and more.
  For details, see the [MP-Opt-Model 4.1 release notes][12] *(also in
  Appendix C in the [MP-Opt-Model User's Manual][13])*.
- [MOST][14] 1.2 adds calculation of expected temporal locational marginal
  price (TLMP), includes transitions into first period in ramping
  reserves, and more. For details, see the [MOST 1.2 release notes][15]
  *(also in Appendix B in the [MOST User's Manual][16])*.
- New options:
  - New AC power flow solver based on `fsolve()` function, selected by
    setting `'pf.alg'` option to `'FSOLVE'`. *(`fsolve()` is part of the MATLAB
    Optimization Toolbox and included in Octave. Supported via MP-Core
    only.)*
  - New Implicit Z-bus Gauss method power flow for distribution systems,
    selected by setting `pf.alg` option to `'ZG'`. *(Support for PV buses
    is very limited.)*
  - New `'exp.use_legacy_core'` option to bypass MP-Core and force use of
    legacy core code for `runpf()`, `runcpf()`, `runopf()`.
- New functions/methods:
  - MP-Test 8.0b1
    - `t_str_match()` -- test that a string matches expected value
    - `t_file_match()` -- test that the content of a text files matches
      that of another file
  - `have_feature_mp_core` - Determines availability of MP-Core.
  - `run_mp` - Top-level function for running any task (PF, CPF, OPF) with
    the new MP-Core and flexible framework.
  - `run_pf` - Wrapper around `run_mp` for running PF.
  - `run_cpf` - Wrapper around `run_mp` for running CPF.
  - `run_opf` - Wrapper around `run_mp` for running OPF.


#### New Case Files:

- Two new European case files. *Thanks to Florin Capitanescu.*
  - `case60nordic` -- 60-bus Nordic case
  - `case8387pegase` -- 8387-bus PEGASE case


#### New Documentation:

- [MATPOWER Developer's Manual][4] -- describes the architecture of the
  new MP-Core and MATPOWER flexible framework
- [*MATPOWER Technical Note 5*][5] "MP-Element: A Unified MATPOWER
  Element Model, with Corresponding Functions and Derivatives"


#### Other Improvements:

- MP-Opt-Model object is now used for power flow as well as OPF and is
  added as `om` field to power flow `results` struct.
- New MATPOWER Docker image (named [`matpower/matpower`][17]) is
  based on the official GNU Octave image ([`gnuoctave/octave`][18]) and
  is available for multiple MATPOWER and Octave versions.
- Each power flow is initialized with the solved voltages of the previous
  one when changing PV buses to PQ during Q limit enforcement.
  *Thanks to Tostado-Véliz, Kamel, Jurado.*
- The `makePTDF()` function can now optionally use a different slack
  distribution for each bus by specifying the `slack` input as a matrix.
  *Thanks to Jon Martinez Corral.*
- Automatically reduce order of polynomial generator costs with higher
  order coefficients equal to zero. Allows the DC OPF to solve cases,
  e.g. with cubic costs where the 3rd order term is 0, such as cases
  exported by PowerWorld. *Thanks to Rajesh Mookerjee.*


#### Bugs Fixed:

- Fatal bug in `int2ext()` when called with `mpopt` and an `int2ext` user
  callback function.
- Bug affecting AC OPF initialization for cases with piecewise linear
  costs and cartesian voltage representations. Fix results in improved
  convergence for affected cases.
- Fix generator voltage set points in `case9target` to match `case9`.
- A vector-valued `label` passed to `apply_changes()` now throws a
  useful error.
- Fix bug in previously undocumented feature of `makePTDF()` where
  the `slack` input is a matrix. *Thanks to Jon Martinez Corral.*


#### Incompatible Changes:

- Removed unused `mpopt` argument from `opf_gen_cost_fcn()` inputs.
- Remove deprecated functions:
  - `d2AIbr_dV2()` -- use `dA2br_dV2()` instead.
  - `d2ASbr_dV2()` -- use `dA2br_dV2()` instead.
  - Deprecated methods of `@opf_model`:
    - `add_constraints()` -- use the corresponding one of the following
      methods instead: `add_lin_constraint()`, `add_nln_constraint()`, or
      `init_indexed_name()`.
    - `add_costs()` -- use the corresponding one of the following methods
      instead: `add_quad_cost()`,  `add_nln_cost()`, `add_legacy_cost()`,
      or `init_indexed_name()`.
    - `add_vars()` -- use the corresponding one of the following methods
      instead: `add_var()`, or `init_indexed_name()`.
    - `build_cost_params()` -- no longer needed, incorporated into
      `params_legacy_cost()`.
    - `get_cost_params()` -- use `params_legacy_cost()` instead.
    - `getv()` -- use `params_var()` instead.
    - `linear_constraints()` -- use `params_lin_constraint()` instead.
- Remove deprecated option:
    - `'opf.init_from_mpc'` -- use `'opf.start'` instead.


[1]: https://github.com/MATPOWER/matpower/blob/master/CHANGES.md
[2]: https://github.com/MATPOWER/matpower/blob/master/docs/MATPOWER-manual.pdf
[3]: https://github.com/MATPOWER/mp-element
[4]: https://matpower.org/documentation/dev-manual/
[5]: https://matpower.org/docs/TN5-MP-Element.pdf
[6]: https://github.com/MATPOWER/mptest
[7]: https://github.com/MATPOWER/mptest/blob/master/docs/relnotes/MP-Test-Release-Notes-8.0.md
[8]: https://github.com/MATPOWER/mips
[9]: https://github.com/MATPOWER/mips/blob/master/docs/relnotes/MIPS-Release-Notes-1.5.md
[10]: https://matpower.org/docs/MIPS-manual-1.5.pdf
[11]: https://github.com/MATPOWER/mp-opt-model
[12]: https://github.com/MATPOWER/mp-opt-model/blob/master/docs/relnotes/MP-Opt-Model-Release-Notes-4.1.md
[13]: https://matpower.org/docs/MP-Opt-Model-manual-4.1.pdf
[14]: https://github.com/MATPOWER/most
[15]: https://github.com/MATPOWER/most/blob/master/docs/relnotes/MOST-Release-Notes-1.2.md
[16]: https://matpower.org/docs/MOST-manual-1.2.pdf
[17]: https://hub.docker.com/r/matpower/matpower
[18]: https://hub.docker.com/r/gnuoctave/octave
