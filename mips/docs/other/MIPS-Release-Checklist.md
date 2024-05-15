MIPS Release Checklist
======================


Pre-release
-----------
- Check [MIPS issue tracker](https://github.com/MATPOWER/mips/issues)
  and `untracked/MIPS-To-Do-List.md` for to do items.
- Create & checkout new `prep-for-release` branch from latest `master`.
- Release notes:
  - Make sure Release History in Appendix C of `MIPS-manual.tex` is
    up-to-date.
  - Create `docs/relnotes/MIPS-Release-Notes-#.#.md` document from
    Appendix C of `MIPS-manual.tex`.
- Update date in Copyright line in:
  - `LICENSE`
  - `docs/sphinx/source/conf.py`.
- Update version number and date in:
  - `mipsver.m`
  - `docs/sphinx/source/conf.py`
  - `lib/Contents.m`
  - `docs/relnotes/MIPS-Release-Notes-#.#.md`
  - `docs/src/MIPS-manual/MIPS-manual.tex`
    - title page
    - copyright (front page and LICENSE text)
    - Appendix C Release History
    - `\mipsver` (update `\mpver`, `\mptestver`, `\mpomver`, `\mostver` too)
  - Sphinx docs
    - `mp-docs-shared/preamble.tex.txt` - \mipsver
    - `mp-docs-shared/prolog.rst.txt` - in URL in raw-html for |MIPSman|
- In `README.md` and `docs/src/MIPS-manual/MIPS-manual.tex`
  - update output of:
    - `test_mips` in Section 2.2
      - (may require `rmpath('/usr/local/pardiso/current')` first)
    - `mips_example2` in Section 3.2 (has version and date)
  - check for any highlighting `\hl`
- Create new DOI for this version of the User's Manual
  - Go to https://doi.org/10.5281/zenodo.3236506
    - Click "New Version" to reserve new DOI for new version
  - Make updates for current version specific citations:
    - version number (3 places)
    - year
    - latest version DOI, current is: 10.5281/zenodo.11176870
      - (update here each time)
    ... in the following places ...
    - CITATION file
    - Citing ... section of README.md
    - Citing ... section of User's Manual
    - Citing ... section of website (not currently here)
  - Make updates for non-version specific citations:
    - search everywhere for 10.5281/zenodo.3236506 and update year
      - User's Manual
      - search citations in all other projects being updated simultaneously
        - MATPOWER-manual.tex
        - MP-Opt-Model-manual.tex
- Create `MIPS-manual.pdf` from `MIPS-manual.tex` and move to `docs`.
- Add release notice with date and version in `CHANGES.md`.
- Commit all changes to `prep-for-release`.
- Push `prep-for-release` to GitHub.
- Make sure CI checks are ok.
- Make copy of `docs/MIPS-manual.pdf` named `MIPS-manual-x.x.pdf`
  - copy to `docs` directory of `matpower.org-static` git repo
    - update `MIPS-manual.pdf` symlink on `https://matpower.org/docs/` to point
      to new `MIPS-manual-x.x.pdf` (replaces existing current version)
      - `cd dev/projects/matpower.org-static/docs`
      - `rm MIPS-manual.pdf`
      - `ln -s ./MIPS-manual-x.x.pdf MIPS-manual.pdf`
    - commit & push, then pull to matpower.org
  - upload `MIPS-manual-x.x.pdf` to Zenodo and finish entry for "New Version"
    - update:
      - Publication date
      - Version
      - Identifiers:
        - version number in "identical to"
- Add link on `https://matpower.org/doc/manuals/` page


Release
-------
- Merge latest `prep-for-release` into `master`.
- Tag with version number, e.g. `1.5`.
- Push `master` to GitHub.
- Publish new release on GitHub: https://github.com/MATPOWER/mips/releases/new
  - use (possibly shortened) contents of `docs/relnotes/MIPS-Release-Notes-#.#.md`


Post-release
------------
- Merge latest `master` into `release`.
- Push `release` to GitHub.
- In manual
  - update version to dev version
  - comment out date line so it uses current date
