MP-Test Release Checklist
=========================


Pre-release
-----------
- Check [MP-Test issue tracker](https://github.com/MATPOWER/mptest/issues)
  for to do items.
- Create & checkout new `prep-for-release` branch from latest `master`.
- Release notes:
  - Create `docs/relnotes/MP-Test-Release-Notes-#.#.md` document
- Update date in Copyright line in:
  - `LICENSE`
  - `docs/sphinx/source/conf.py`.
- Update version number and date in:
  - `mptestver.m`
  - `docs/sphinx/source/conf.py`.
- Pull in latest `mp-docs-shared` and build and check Sphinx docs.
- Update output in `README.md`, as necessary, of:
  - `test_mptest`
  - `mptest_ex1`
  - `test_everything_ex1`
- Add release notice with date and version in `CHANGES.md`.
- Add/update (including release date) docs/relnotes/MP-Test-Release-Notes-x.x.md
- Commit all changes to `prep-for-release`.
- Push `prep-for-release` to GitHub.
- Make sure CI checks are ok.

Release
-------
- Merge latest `prep-for-release` into `master`.
- Tag with version number, e.g. `8.0`.
- Push `master` to GitHub.
- Publish new release on GitHub: https://github.com/MATPOWER/mptest/releases/new

Post-release
------------
- Merge latest `master` into `release`.
- Push `release` to GitHub.
