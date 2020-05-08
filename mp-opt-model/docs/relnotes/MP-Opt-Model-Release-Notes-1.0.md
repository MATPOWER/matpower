What's New in MP-Opt-Model 1.0
------------------------------

#### Released May 8, 2020

Below is a summary of the changes since version 0.8 of MP-Opt-Model. See
the [`CHANGES.md`][1] file for all the gory details. For release notes
for previous versions, see Appendix C of the [MP-Opt-Model User's
Manual][2].

#### New Documentation:
  - Add the [MP-Opt-Model User's Manual][2] with LaTeX source code
    included in `docs/src`.

#### Other Improvements:
  - Refactor `opt_model` class to inherit from new abstract base class
    `mp_idx_manager` which can be used to manage the indexing of other
    sets of parameters, etc. in other contexts.


[1]: https://github.com/MATPOWER/mp-opt-model/blob/master/CHANGES.md
[2]: https://github.com/MATPOWER/mp-opt-model/blob/master/docs/MP-Opt-Model-manual.pdf
