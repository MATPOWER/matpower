MATPOWER Release Checklist
==========================

Pre-release
-----------

- Prepare release landing page on website.
  - Create new draft download to get download IDs.
  - Use download ID to create new download button widget.
- Prepare release e-mail(s).
- Check value of `min_ver.Octave` and `min_ver.MATLAB` in `install_matpower.m`
- Check [MATPOWER issue tracker](https://github.com/MATPOWER/matpower/issues)
  and `untracked/MATPOWER-To-Do-List.md` for to do items.
- Prepare `mp-docs-shared`.
  - Update date in copyright line in `mp-docs-shared/LICENSE`.
  - Update all version numbers.
    - in `preamble.tex.txt`, definitions of:
      - `\mptestver`
      - `\mipsver`
      - `\mpomver`
      - `\mostver`
      - `\mpver`
    - in `prolog.rst.txt`, in URLs in definitions of:
      - `|MUM|`
      - `|MIPSman|`
      - `|MPOMman|`
      - `|MOSTman|`
- Release MP-Test.
- Release MIPS.
- Release MP-Opt-Model.
- Release MOST.
- Release MATPOWER Extras.
- Release BPMPD (if necessary, check that it works).
- Release MINOPF (if necessary, check that it works).
- Release TSPOPF (if necessary, check that it works).
- Create & checkout new `prep-for-release` branch from latest `master`.
- Get released subrepos (use `--branch=release` if different from `master`):
  - `git subrepo pull --branch=master docs/sphinx/source/mp-docs-shared`
  - `git subrepo pull --branch=master mptest`
  - `git subrepo pull --branch=master mips`
  - `git subrepo pull --branch=master mp-opt-model`
  - `git subrepo pull --branch=master most`
- Release notes:
  - Make sure Release History in Appendix H of `MATPOWER-manual.tex` is
    up-to-date.
  - Create `docs/relnotes/MATPOWER-Release-Notes-#.#.md` document from
    Appendix H of `MATPOWER-manual.tex`.
  - Create `docs/relnotes/MATPOWER-Announce-#.#.md` document from
    `docs/relnotes/MATPOWER-Release-Notes-#.#.md`.
- Update version number and date in:
  - `mpver.m`
  - `lib/Contents.m`
  - `docs/relnotes/MATPOWER-Release-Notes-#.#.md`
  - `docs/relnotes/MATPOWER-Announce-#.#.md`
  - `docs/src/MATPOWER-manual/MATPOWER-manual.tex`
    - title page
    - copyright (front page and LICENSE text)
    - Appendix H Release History
    - `\mpver` (update `\mptestver`, `\mpomver`, `\mipsver`, `\mostver`, `\sdpopfver` and `\syngridver` too)
    - do search for all possible old versions 6.1-dev, 7.0b1
  - Copyright line in `LICENSE`.
- In `docs/src/MATPOWER-manual/MATPOWER-manual.tex`
  - update output of:
    - `test_matpower` (Section 2.3, Step 3), use `test_matpower_no_options`
    - `help runopf` (Section 2.4)
    - `runcpf` (Section 5.6.2, Example)  (run `cpf_example`)
  - check for any highlighting `\hl`
  - check for pagination anomalies (use \clearpage)
- Create new DOI for this version of the User's Manual
  - Go to https://doi.org/10.5281/zenodo.3236519
    - Click "New Version" to reserve new DOI for new version
  - Make updates for current version specific citations:
    - version number (3 places)
    - year
    - latest version DOI, current is: 10.5281/zenodo.11212313
      - (update here each time)
    ... in the following places ...
    - CITATION file
    - Citing ... section of README.md
    - Citing ... section of User's Manual
    - Citing ... section of website (https://matpower.org/citing/)
  - Make updates for non-version specific citations:
    - search everywhere for 10.5281/zenodo.3236519 and update year
      - MATPOWER User's Manual
      - search citations in all other projects being updated simultaneously
        - e.g. MOST, SynGrid, MP-Sim, sopf1 manuals, TN2, TN3, TN4
- Create new DOI for this version of the software
  - Go to https://doi.org/10.5281/zenodo.3236535
    - Click "New Version" to reserve new DOI for new version
  - Make updates for current version specific citations:
    - version number (2 places)
    - year
    - latest version DOI, current is: 10.5281/zenodo.11212330
      - (update here each time)
    ... in the following places ...
    - CITATION file
    - Citing ... section of README.md
    - Citing ... section of User's Manual
    - Citing ... section of website (https://matpower.org/citing/)
  - Make updates for non-version specific citations:
    - search everywhere for 10.5281/zenodo.3236535 and update year
      - User's Manual
      - search citations in all other projects being updated simultaneously
        - e.g. MIPS, MOST, SynGrid, MP-Sim, sopf1 manuals, TN2, TN3, TN4
- Copy latest `MIPS-manual.aux`, `MP-Opt-Model-manual.aux`, `MOST-manual.aux`
  to `docs/src/MATPOWER-manual` for `\externaldocument`
- Create `MATPOWER-manual.pdf` from `MATPOWER-manual.tex` and move to `docs`.
- Update Docker files
  - follow instructions in `docker/Docker-Build-Notes.md` for new Octave
    images, if applicable
  - update value listed for _current release_ of MATPOWER in **Versions** section
    of `docker/MATPOWER-Docker.md`, and add new rows to table.
  - add lines in `build_matpower_images.sh` for new MATPOWER version
  - build and test a new "local" image with MATPOWER Extras installed
- Add release notice with date and version in `CHANGES.md`.
- Commit all changes to `prep-for-release`.
- Push `prep-for-release` to GitHub.
- Make sure CI checks are ok.
- Make copy of `docs/MATPOWER-manual.pdf` named `MATPOWER-manual-x.x.pdf`
  - copy to `docs` directory of `matpower.org-static` git repo
    - update `MATPOWER-manual.pdf` symlink on `https://matpower.org/docs/` to point
      to new `MATPOWER-manual-x.x.pdf` (replaces existing current version)
      - `cd dev/projects/matpower.org-static/docs`
      - `rm MATPOWER-manual.pdf`
      - `ln -s ./MATPOWER-manual-x.x.pdf MATPOWER-manual.pdf`
    - commit & push, then pull to matpower.org
  - upload `MATPOWER-manual-x.x.pdf` to Zenodo and finish entry for "New Version"
    - update:
      - Publication date
      - Version
      - Identifiers:
        - version number in "identical to"
        - version specific DOI for "Documents"
  - add link on `https://matpower.org/doc/manuals/` page


Release
-------
- Merge latest `prep-for-release` into `master`.
- Tag with version number, e.g. `8.0`.
- Push `master` to GitHub.
- create archive with MATPOWER-Extras
  - `cd matpower/untracked/`
  - `ship_it.sh 7.0b1`
  - (obsolete) copy resulting `matpower#.#.zip` file to
    `https://matpower.org/downloads/#.#/`
- upload to Zenodo and finish entry for "New Version"
  - update:
    - Publication date
    - Version
    - Identifiers:
      - version specific DOI for "documents this upload" (MATPOWER manual)
      - version specific DOI for "documents this upload" (MOST manual)
      - version specific DOI for "documents this upload" (MIPS manual)
      - version specific DOI for "documents this upload" (MP-Opt-Model manual)
- Publish new release on GitHub: https://github.com/MATPOWER/most/releases/new
  - use contents of `docs/relnotes/MATPOWER-Announce-#.#.md`
  - add download badge at top:
    - ![7.0 Downloads](https://img.shields.io/github/downloads/MATPOWER/matpower/7.0/total.svg)
  - add DOIs at bottom:
    -   General DOI:
        [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3236535.svg)](https://doi.org/10.5281/zenodo.3236535)

        Version Specific DOI:
        [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3251119.svg)](https://doi.org/10.5281/zenodo.3251119)
- (obsolete)
  - Add version number to
    - `https://matpower.org/sw-download-config.php`
    - `https://matpower.org/index.html`
- Update matpower.org website
  - create new "download" pointing to GitHub release
    - use link to GitHub release, don't upload the file
      - this gives a size of ~640 bytes, which has to be corrected manually
        in the database (since it displays in some tables), use phpMyAdmin or
        other db editor, find the row in the `wp_postmeta` table with
        `post_id` equal to the download id for the file, and
        `meta_key` equal to `_filesize`, set the value equal to the file
        size in bytes
        (note that the post_id may be 1 more than the download ID)
    - create new Custom HTML widget for new Download button
    - use this widget in release landing page
    - update download IDs in release landing page
    - remove Current tag from previous release and add it to new release
      - not sure I'm actually using tags
    - remove "Current" category from previous release and set it to "MATPOWER"
      - add "Current" category to new release
  - update "Downloads" sidebar widget as necessary (e.g. no pre-release avaiable)
  - create new blog post announcing release (copied from docs/relnotes/MATPOWER-Announce-#.#.md)
    - use "Releases" category
  - update download counter file at rdzman@matpower.org:/home4/rdzman/matpower.org-downloads
    - add a column to the header line at the top of matpower-download-log.txt
  - build latest MATPOWER Sphinx documentation
    - commit to matpower-docs
    - push
    - pull to matpower.org
- MATLAB Central File Exchange
    - https://www.mathworks.com/matlabcentral/fileexchange/72085-matpower/
    (this is currently just a link)
- Build new Docker images according to instructions in
  `docker/Docker-Build-Notes.md` for new MATPOWER images


Post-release
------------
- Merge latest `master` into `release`.
- Push `release` to GitHub.
- In manual
  - update version to dev version
  - comment out date line so it uses current date
- new dev version number
- update web site
- `cd /usr/local/matpower/dist`
- `git clone --branch 7.1 git@github.com:MATPOWER/matpower.git matpower-version7_1`
- `cd /usr/local/matpower/dist/matpower-version7_1`
- `rm -rf .git`
- `git clone --branch 7.1 git@github.com:MATPOWER/matpower-extras.git extras`
- `cd extras`
- `rm -rf .git`


To create the Online Function Reference
---------------------------------------

- unzip in ~/dev/temp
- delete `data` subdirectory
- in Matlab:
```matlab
    cd ~/dev/temp
    v = '7.1';
    m2html('mfiles',['matpower' v],'htmldir',['htdocs/docs/ref' v],'global','on','recursive','on','index',['menu' v],'template','matpower', 'graph', 'off')
```
- then:
```
    mv ~/dev/temp/htdocs/docs/ref7.1/menu7.1.html ~/dev/projects/matpower.org-static/docs/ref/
    mv ~/dev/temp/htdocs/docs/ref7.1/matpower7.1 ~/dev/projects/matpower.org-static/docs/ref/
```
  - update index.html symbolic link to latest version
```
    cd ~/dev/projects/matpower.org-static/docs/ref/
    rm index.html
    ln -s ./menu7.1.html index.html
```
  - delete beta versions, e.g. 7.0b1, then symlink them
    so old direct links point to the 7.1 release
- commit everything
  - do change of `index.html` as a separate 2nd commit (and pulled separately too)
    - otherwise pull can fail if it tries to update `index.html` before creating
      the target file





---------

Old
---

    commit all changes
    push all changes
    build the dist archive:
        git clone git@bitbucket.org:rdzman/matpower.git
        cd matpower
        git checkout develop
        mv dist matpower6.0
        zip -9 -r matpower6.0.zip matpower6.0
    put it in new download dir with license file
    unzip in ~/dev/temp and create HTML docs:
        in Matlab:
            cd ~/dev/temp
            v = '6.0';
            m2html('mfiles',['matpower' v],'htmldir',['htdocs/docs/ref' v],'global','on','recursive','on','index',['menu' v],'template','matpower', 'graph', 'off')
        then:
            mv ~/dev/temp/htdocs/docs/ref6.0/menu6.0.html ~/dev/projects/matpower/htdocs/docs/ref/
            mv ~/dev/temp/htdocs/docs/ref6.0/matpower6.0 ~/dev/projects/matpower/htdocs/docs/ref/
            update index.html symbolic link to latest version
                cd ~/dev/projects/matpower/htdocs/docs/ref/
                rm index.html
                ln -s ./menu6.0.html index.html
            delete beta versions, e.g. 6.0b1 and 6.0b2, then symlink them
                so old direct links point to the 6.0 release
    commit everything and tag it, merge to master branch
    check out into /usr/local/matpower/dist
        git clone git@bitbucket.org:rdzman/matpower.git
        cd matpower
        git checkout 6.0
        mv dist ../matpower-version6_0
        cd ..
        rm -rf /usr/local/matpower/dist/matpower
    create /usr/local/matpower/v6.0

---------

    old CVS/SVN-based commands:
        build the dist archive:
            cvs ex -D now matpower
            svn export http://shrike.calsnet.cornell.edu:9880/psercsvn/matpower/trunk matpower
            mv matpower/dist matpower5.0
            svn export http://shrike.calsnet.cornell.edu:9880/psercsvn/matpower/trunk/dist matpower5.0
        commit everything and tag it, merge to master branch
            svn copy http://shrike.calsnet.cornell.edu:9880/psercsvn/matpower/trunk \
            http://shrike.calsnet.cornell.edu:9880/psercsvn/matpower/tags/version5_1 \
            -m "Release 5.0"
        check out into /usr/local/matpower/dist
            svn co http://shrike.calsnet.cornell.edu:9880/psercsvn/matpower/tags/version5_1/dist matpower-version5_1
