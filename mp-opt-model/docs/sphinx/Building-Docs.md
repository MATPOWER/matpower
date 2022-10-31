How to Build the Documentation
==============================

The documentation for this project is based on [Sphinx][1], a [Python][2]-based
documentation processing package. The documentation source is written in the
[reStructuredText][3] (rST) format and compiled to HTML and/or PDF formats.

As part of [MATPOWER][4], the documentation for this project is normally
built as a component of the overall MATPOWER documentation. However, the HTML
version can be build separately by following the steps below.

1.  Install [Python][2].

2.  Create and activate a Python [virtual environment][5]. *(optional, but recommended)*

    Unix/macOS
    ```
    python3 -m venv sphinx-env
    source sphinx-env/bin/activate
    ```

    Windows
    ```
    python3 -m venv sphinx-env
    .\sphinx-env\Scripts\activate
    ```

    Then confirm you're in the virtual environment:

    Unix/macOS
    ```
    which python
    ```

    Windows
    ```
    where python
    ```

    To leave the virtual environment:
    ```
    deactivate
    ```

3.  Install [Sphinx][1], [sphinxcontrib-matlabdomain][6], [sphinx-tabs][7] and
    [sphinx-rtd-theme][8].
    ```
    pip install -U sphinx
    pip install -U sphinxcontrib-matlabdomain
    pip install -U sphinx-tabs
    pip install -U sphinx-rtd-theme
    ```

4.  Include [mp-docs-shared][9], a collection of files shared by all of the
    [Sphinx][1] documention for [MATPOWER][4], in the `source` directory.
    ```
    cd docs/sphinx/source
    git clone https://github.com/MATPOWER/mp-docs-shared.git
    ```

5.  In the `docs/sphinx` directory, type:
    ```
    make html
    ```

The resulting HTML documentation can be found in `docs/sphinx/build/html` and
accessed by opening `docs/sphinx/build/html/index.html` in a web browser.


----
[1]: https://www.sphinx-doc.org
[2]: https://python.org
[3]: https://www.writethedocs.org/guide/writing/reStructuredText/
[4]: https://matpower.org
[5]: https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment
[6]: https://pypi.org/project/sphinxcontrib-matlabdomain/
[7]: https://pypi.org/project/sphinx-tabs/
[8]: https://pypi.org/project/sphinx-rtd-theme/
[9]: https://github.com/MATPOWER/mp-docs-shared/
