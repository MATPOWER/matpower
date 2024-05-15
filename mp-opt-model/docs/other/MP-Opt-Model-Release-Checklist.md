MP-Opt-Model Release Checklist
==============================


Pre-release
-----------
- Check [MP-Opt-Model issue tracker](https://github.com/MATPOWER/mp-opt-model/issues)
  and `untracked/MP-Opt-Model-To-Do-List.md` for to do items.
- Create & checkout new `prep-for-release` branch from latest `master`.
- Release notes:
  - Make sure Release History in Appendix C of `MP-Opt-Model-manual.tex` is
    up-to-date.
  - Create `docs/relnotes/MP-Opt-Model-Release-Notes-#.#.md` document from
    Appendix C of `MP-Opt-Model-manual.tex`.
- Update date in Copyright line in:
  - `LICENSE`
  - `docs/sphinx/source/conf.py`.
- Update version number and date in:
  - `mpomver.m`
  - `docs/sphinx/source/conf.py`
  - `lib/Contents.m`
  - `docs/relnotes/MP-Opt-Model-Release-Notes-#.#.md`
  - `docs/src/MP-Opt-Model-manual/MP-Opt-Model-manual.tex`
    - title page
    - copyright (front page and LICENSE text)
    - Appendix C Release History
    - `\mpomver` (update `\mpver`, `\mptestver`, `\mostver`, `mipsver` too)
  - Sphinx docs
    - `mp-docs-shared/preamble.tex.txt` - \mpomver
    - `mp-docs-shared/prolog.rst.txt` - in URL in raw-html for |MPOMman|
- In `README.md` and `docs/src/MP-Opt-Model-manual/MP-Opt-Model-manual.tex`
  - update output of:
    - `test_mp_opt_model` in Section 2.2
    - `nlps_master_ex2` in Section 4.3.2 (has MIPS version and date)
    - `pne_ex1` in Section 4.5.8 (has MP-Opt-Model version and date)
  - check for any highlighting `\hl`
- Create new DOI for this version of the User's Manual
  - Go to https://doi.org/10.5281/zenodo.3818002
    - Click "New Version" to reserve new DOI for new version
  - Make updates for current version specific citations:
    - version number (3 places)
    - year
    - latest version DOI, current is: 10.5281/zenodo.11177079
      - (update here each time)
    ... in the following places ...
    - CITATION file
    - Citing ... section of README.md
    - Citing ... section of User's Manual
    - Citing ... section of website (not currently here)
  - Make updates for non-version specific citations:
    - search everywhere for 10.5281/zenodo.3818002 and update year
      - User's Manual
      - search citations in all other projects being updated simultaneously
        - MATPOWER User's Manual
- Copy latest `MIPS-manual.aux` to `docs/src/MP-Opt-Model-manual` for
  `\externaldocument`
- Create `MP-Opt-Model-manual.pdf` from `MP-Opt-Model-manual.tex` and move
  to `docs`.
- Add release notice with date and version in `CHANGES.md`.
- Commit all changes to `prep-for-release`.
- Push `prep-for-release` to GitHub.
- Make sure CI checks are ok.
- Make copy of `docs/MP-Opt-Model-manual.pdf` named `MP-Opt-Model-manual-x.x.pdf`
  - copy to `docs` directory of `matpower.org-static` git repo
    - update `MP-Opt-Model-manual.pdf` symlink on `https://matpower.org/docs/` to point
      to new `MP-Opt-Model-manual-x.x.pdf` (replaces existing current version)
      - `cd dev/projects/matpower.org-static/docs`
      - `rm MP-Opt-Model-manual.pdf`
      - `ln -s ./MP-Opt-Model-manual-x.x.pdf MP-Opt-Model-manual.pdf`
    - commit & push, then pull to matpower.org
  - upload `MP-Opt-Model-manual-x.x.pdf` to Zenodo and finish entry for "New Version"
    - update:
      - Publication date
      - Version
      - Identifiers:
        - version number in "identical to"
- Add link on `https://matpower.org/doc/manuals/` page


Release
-------
- Merge latest `prep-for-release` into `master`.
- Tag with version number, e.g. `4.2`.
- Push `master` to GitHub.
- Publish new release on GitHub: https://github.com/MATPOWER/mp-opt-model/releases/new
  - use (possibly shortened) contents of `docs/relnotes/MP-Opt-Model-Release-Notes-#.#.md`
- MATLAB Central File Exchange
    - https://www.mathworks.com/matlabcentral/fileexchange/77909-mp-opt-model
    (this is currently a GitHub linked entry)


Post-release
------------
- Merge latest `master` into `release`.
- Push `release` to GitHub.
- In manual
  - update version to dev version
  - comment out date line so it uses current date
