What's New in MATPOWER 7.0
--------------------------

#### Released May 17, 2024

Below is a summary of the changes since version 8.0b1 of MATPOWER. See the
[`CHANGES.md`][1] file for all the gory details. For release notes for
previous versions, see Appendix H of the [MATPOWER User's Manual][2].

#### New Features:
- [MP-Opt-Model][11] 4.2 adds support for:
  - Memory savings via optional use of transpose of linear constraint matrix.
  - Display of solved model solution.
  - MATLAB Optimization Toolbox 9.5 and later (MATLAB R2023a and later).

  For details, see the release notes for MP-Opt-Model [4.2][12a] *(also in
  Appendix C in the [MP-Opt-Model User's Manual][13])*.
- [MOST][14] 1.3 improves speed and memory usage on certain problems. For
  details, see the release notes for MOST [1.3][15a] *(also
  in Appendix B in the [MOST User's Manual][16])*.
- [SimulinkMATPOWER][16a], included in [MATPOWER Extras][16b], enables the use
  of MATPOWER from MATLAB's [Simulink®][16c] environment.
  *Thanks to Lukas Ortmann.*


#### New Case Files:

- Two new Swedish distribution system cases. *Thanks to Gabriel Malmer.*
  - `case533mt_lo` -- 533-bus case, DSO Kraftringen, *low* net-load scenario
  - `case533mt_hi` -- 533-bus case, DSO Kraftringen, *high* net-load scenario


#### New Documentation:

- [MATPOWER Documentation site][16d] -- new documentation website with MATPOWER
  manuals, reference docs, and how-to guides
- [MATPOWER Reference Manual][4a] -- reference for MATPOWER functions and
  classes, especially those in the new MP-Core and MATPOWER flexible framework


#### Other Improvements:
- New option to `makeLODF()` to replace with `NaN` any columns corresponding
  to "bridge" branches which, if removed, result in islands.
  *Thanks to Liangyu Zhang.*
- Update versions of included packages:
  - MIPS 1.5.1
  - MP-Opt-Model 4.2
  - MOST 1.3
  - MP-Test 8.0

#### Bugs Fixed:
- Fix bug #210 where radial power flow methods resulted in numerical error
  with multiple generators at a bus. *Thanks to Mirko Todorovski.*
- Fix bug #223, fatal error in `save2psse()` for cases with dispatchable loads.
  *Thanks to Ali Jahanbani.*
- Fix crash in `radial_pf()` under MATLAB R2016a and earlier.
- Delete extra (duplicate) row in `gencost` in `case94pi.m`.


[1]: https://github.com/MATPOWER/matpower/blob/master/CHANGES.md
[2]: https://github.com/MATPOWER/matpower/blob/master/docs/MATPOWER-manual.pdf
[4a]: https://matpower.org/doc/ref-manual/
[11]: https://github.com/MATPOWER/mp-opt-model
[12a]: https://github.com/MATPOWER/mp-opt-model/blob/master/docs/relnotes/MP-Opt-Model-Release-Notes-4.2.md
[13]: https://matpower.org/docs/MP-Opt-Model-manual.pdf
[14]: https://github.com/MATPOWER/most
[15a]: https://github.com/MATPOWER/most/blob/master/docs/relnotes/MOST-Release-Notes-1.3.md
[16]: https://matpower.org/docs/MOST-manual.pdf
[16a]: https://github.com/MATPOWER/mx-simulink_matpower
[16b]: https://github.com/MATPOWER/matpower-extras
[16c]: https://www.mathworks.com/products/simulink.html
[16d]: https://matpower.org/doc/
