What's New in MATPOWER 7.0
--------------------------

#### Released Jun 20, 2019

Below is a summary of the changes since version 7.0b1 of MATPOWER. See the
[`CHANGES.md`][1] file for all the gory details. For release notes for
previous versions, see Appendix H of the [MATPOWER User's Manual][2].

#### New Features:
- Three new variants of the standard Newton AC power flow, for a total
  of four, including both nodal power and current balance constraints
  and both polar and cartesian representations of voltage. See the new
  `pf.current_balance` and `pf.v_cartesian` options.
  *Thanks to Baljinnyam Sereeter.*
- Docker image tagged [`matpower/matpower-desktop`][10], providing a
  pre-packaged Ubuntu desktop environment with Octave, MATPOWER, and
  the [MATPOWER Extras][11] all pre-installed. See the
  [`docker/MATPOWER-Docker.md`][12] page for more details.
  *Thanks to Richard Lincoln.*
- New options:
  - `pf.current_balance` and `pf.v_cartesian` control formulation used
    for Newton AC power flow. *Thanks to Baljinnyam Sereeter.*

#### New Website
- MATPOWER has a new website at https://matpower.org. Please discontinue
  use of the old http://www.pserc.cornell.edu/matpower/ address.

#### Other Improvements:
- Update versions of included packages:
  - MIPS 1.3.1.
  - MOST 1.0.2.
  - MP-Test 7.0.
- Improve performance of `makeYbus()`. *Thanks to Binbin Chen.*
- Add support for IPOPT solver under Octave, including in the Travis-CI
  testing.
  *Thanks to Richard Lincoln.*
- Add support for YALMIP, SeDuMi and SDPT3 to be recognized under Octave.
  *Thanks to Richard Lincoln.*
- Add `'clear_cache'` options to `have_fcn()` (see #65) to facilitate
  re-checking for optional functionality after changes to the MATLAB/Octave
  path.

#### Bugs Fixed:
- Fix #53 where certain OPF cases (e.g. `case33bw`) resulted in a fatal
  error under versions of MATLAB prior to R2016b (v9.1).
  *Thanks to Jane Cheung.*
- Fix #56 where `save2psse` was missing entries for two transformer
  columns, namely, `VMA1` and `VMI1`.
  *Thanks to Amin Gholami.*
- Fix #57 where `save2psse` always used 1 for the `CKT` number, even for
  parallel branches or transformers.
  *Thanks to Amin Gholami.*
- Fix bug in `have_fcn()` where it would incorrectly mark Gurobi as
  available even if it had an expired license or failed for some other
  reason.
- Fix #60 by adding missing generator at slack bus in RTE cases. Affects
  the following cases:
  - `case1888rte`
  - `case1951rte`
  - `case2848rte`
  - `case2868rte`
  - `case6468rte`
  - `case6470rte`
  - `case6495rte`
  - `case6515rte`
  *Thanks to Jean Maeght.*

#### Incompatible Changes:
- Eliminate unnecessary reordering of on-line generators (sorted by
  increasing bus index) from `ext2int()`. The order is now left
  unmodified by `ext2int()`. This change should only affect user code
  that explicitly depends on the order of generators with internal
  numbering (hopefully quite rare).


[1]: https://github.com/MATPOWER/matpower/blob/master/CHANGES.md
[2]: https://github.com/MATPOWER/matpower/blob/master/docs/MATPOWER-manual.pdf
[3]: https://github.com/GridMod/RTS-GMLC
[4]: https://matpower.org
[5]: https://matpower.org/docs/TN3-More-OPF-Derivatives.pdf
[6]: https://matpower.org/docs/TN4-OPF-Derivatives-Cartesian.pdf
[7]: https://matpower.org/docs/MATPOWER-manual-7.0.pdf
[8]: https://matpower.org/docs/MIPS-manual-1.3.pdf
[9]: https://matpower.org/docs/MOST-manual-1.0.1.pdf
[10]: https://hub.docker.com/r/matpower/matpower-desktop
[11]: https://github.com/MATPOWER/matpower-extras
[12]: https://github.com/MATPOWER/matpower/blob/master/docker/MATPOWER-Docker.md
