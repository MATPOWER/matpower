What's New in MATPOWER 7.0b1
----------------------------

#### Released Oct 31, 2018

Below are some of the highlights of the changes since version 6.0 of
MATPOWER. See the [full release notes][1] and the [`CHANGES.md`][2]
file for more details. For release notes for previous versions, see
Appendix H of the [MATPOWER User's Manual][3].

#### New Features:
- New MATPOWER installer script `install_matpower()`
- User-defined general nonlinear constraints and costs in AC OPF
- PSS/E RAW export
- Cartesian coordinate voltage and current mismatch formulations of AC OPF.
- Three new radial power flow algorithms
- Major update to OPF soft limit functionality
- Many new functions and program options.

#### New Case Files:
- Seven new purely synthetic cases, up to 82,000 buses.
- New [RTS-GMLC][4] case.
- Six new radial distribution system cases.

#### New Documentation:
- Two new Tech Notes, [TN3][5] and [TN4][6]
- LaTeX source code for manuals

#### Other Improvements:
- Update versions of included packages:
  - MIPS 1.3.
  - MOST 1.0.1.
  - MP-Test 7.0b1.
- Continuous integration testing via GitHub and Travis-CI integration.
- Support in core optimization model for:
  - general nonlinear constraints
  - general nonlinear costs
  - quadratic costs
- Refactor OPF code to take advantage of new `opt_model` capabilities
  for nonlinear constraints and quadratic and nonlinear costs.
- Support for polar and cartesian voltages in derivative functions.
- Improved performance (up to 2x speedup) for Newton power flow.
- Handling of generator types, fuel types and bus names.
- Numerous bug fixes.

#### Incompatible Changes:
- Move included MATPOWER case files to new `data` subdirectory.
- Default soft limit behavior relaxes all constraints.
- Minor corrections to data for Polish system cases.
- Add `mpopt` to input args for some OPF callbacks.


[1]: https://github.com/MATPOWER/matpower/blob/master/docs/relnotes/MATPOWER-Release-Notes-7.0.md
[2]: https://github.com/MATPOWER/matpower/blob/master/CHANGES.md
[3]: https://github.com/MATPOWER/matpower/blob/master/docs/MATPOWER-manual.pdf
[3]: https://github.com/GridMod/RTS-GMLC
[4]: http://www.pserc.cornell.edu/matpower/TN3-More-OPF-Derivatives.pdf
[5]: http://www.pserc.cornell.edu/matpower/TN4-OPF-Derivatives-Cartesian.pdf
