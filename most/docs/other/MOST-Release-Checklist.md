MOST Release Checklist
======================


Pre-release
-----------
- Check [MOST issue tracker](https://github.com/MATPOWER/most/issues)
  and `untracked/MOST-To-Do-List.md` for to do items.
- Create & checkout new `prep-for-release` branch from latest `master`.
- Release notes:
  - Make sure Release History in Appendix B of `MOST-manual.tex` is
    up-to-date.
  - Create `docs/relnotes/MOST-Release-Notes-#.#.md` document from
    Appendix B of `MOST-manual.tex`.
- Update date in Copyright line in:
  - `LICENSE`
  - `docs/sphinx/source/conf.py`.
- Update version number and date in:
  - `mostver.m`
  - `docs/sphinx/source/conf.py`
  - `most.m`    (copyright line in printout)
  - `lib/Contents.m`
  - `docs/relnotes/MOST-Release-Notes-#.#.md`
  - `docs/src/MOST-manual/MOST-manual.tex`
    - title page
    - copyright (front page and LICENSE text)
    - Appendix B Release History
    - `\mostver`  (update `\mpver`, `\mptestver`, `\mipsver`, `\mpomver` too)
  - Sphinx docs
    - `mp-docs-shared/preamble.tex.txt` - \mostver
    - `mp-docs-shared/prolog.rst.txt` - in URL in raw-html for |MOSTman|
- In `README.md` and `docs/src/MOST-manual/MOST-manual.tex`
  - update output of `test_most` in Section 2.2
    - may require first doing:
      ```matlab
        have_feature('cplex', 0);
        have_feature('glpk', 0);
        have_feature('intlinprog', 0);
        have_feature('mosek', 0);
        rmpath('/Users/ray/dev/projects/sopf/dist');
      ```
- In `docs/src/MOST-manual/MOST-manual.tex`
  - update output of run of sample case in Section 2.3.2
      ```matlab
        mpc = loadcase('ex_case3b');
        transmat = ex_transmat(12);
        xgd = loadxgendata('ex_xgd_uc', mpc);
        [iwind, mpc, xgd] = addwind('ex_wind_uc', mpc, xgd);
        [iess, mpc, xgd, sd] = addstorage('ex_storage', mpc, xgd);
        contab = ex_contab();
        profiles = getprofiles('ex_load_profile');
        profiles = getprofiles('ex_wind_profile', profiles, iwind);
        mdi = loadmd(mpc, transmat, xgd, sd, contab, profiles);
        mpopt = mpoption('verbose', 1);
        mdo = most(mdi, mpopt);
      ```
  - check for any highlighting `\hl`
- Create new DOI for this version of the User's Manual
  - Go to https://doi.org/10.5281/zenodo.3236531
    - Click "New Version" to reserve new DOI for new version
  - Make updates for current version specific citations:
    - version number (3 places)
    - year
    - latest version DOI, current is: 10.5281/zenodo.11177189
      - (update here each time)
    ... in the following places ...
    - CITATION file
    - Citing ... section of README.md
    - Citing ... section of User's Manual
    - Citing ... section of website (not currently here)
  - Make updates for non-version specific citations:
    - search everywhere for 10.5281/zenodo.3236531 and update year
      - User's Manual
      - search citations in all other projects being updated simultaneously
        - MATPOWER User's Manual
    - search everywhere for 10.5281/zenodo.3236519 and update year (MATPOWER User's Manual)
      - User's Manual
- Copy latest `MATPOWER-manual.aux` to `docs/src/MOST-manual` for
  `\externaldocument`
- Create `MOST-manual.pdf` from `MOST-manual.tex` and move to `docs`.
- Add release notice with date and version in `CHANGES.md`.
- Commit all changes to `prep-for-release`.
- Push `prep-for-release` to GitHub.
- Make sure CI checks are ok.
- Make copy of `docs/MOST-manual.pdf` named `MOST-manual-x.x.pdf`
  - copy to `docs` directory of `matpower.org-static` git repo
    - update `MOST-manual.pdf` symlink on `https://matpower.org/docs/` to point
      to new `MOST-manual-x.x.pdf` (replaces existing current version)
      - `cd dev/projects/matpower.org-static/docs`
      - `rm MOST-manual.pdf`
      - `ln -s ./MOST-manual-x.x.pdf MOST-manual.pdf`
    - commit & push, then pull to matpower.org
  - upload `MOST-manual-x.x.pdf` to Zenodo and finish entry for "New Version"
    - update:
      - Publication date
      - Version
      - Identifiers:
        - version number in "identical to"
- Add link on `https://matpower.org/doc/manuals/` page


Release
-------
- Merge latest `prep-for-release` into `master`.
- Tag with version number, e.g. `1.3`.
- Push `master` to GitHub.
- Publish new release on GitHub: https://github.com/MATPOWER/most/releases/new
  - use (possibly shortened) contents of `docs/relnotes/MOST-Release-Notes-#.#.md`


Post-release
------------
- Merge latest `master` into `release`.
- Push `release` to GitHub.
- In manual
  - update version to dev version
  - comment out date line so it uses current date
