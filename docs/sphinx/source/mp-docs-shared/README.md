mp-docs-shared
==============

### Files shared by multiple [MATPOWER][1] documentation projects.

[MATPOWER][1] documentation projects, based on [Sphinx][2], expect to find
an `mp-docs-shared` directory with the contents of this repo in their
`source` directory. The main [MATPOWER][1] project already includes it as
a [`git subrepo`][3]. However, building the documentation for individual
components separetely (e.g. MP-Test, MIPS, MP-Opt-Model, etc.) requires
that `mp-docs-shared` be manually included in the `source` directory before
building.
```
cd <project>/docs/sphinx/source
git clone https://github.com/MATPOWER/mp-docs-shared.git
```

Each project typically also relies on the use of a symlink for the CSS file.
```
cd <project>/docs/sphinx/source/_static/css
ln -s ../../mp-docs-shared/css/matpower.css matpower.css
```


[1]: https://matpower.org
[2]: https://www.sphinx-doc.org/
[3]: https://github.com/ingydotnet/git-subrepo
