How To Build the Documentation
==============================

1. Install `Python <https://python.org>`_.

2. Create and activate a `virtual environment <https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment>`_ *(optional, but recommended)*

   .. tabs::

      .. group-tab:: Unix/macOS
      
         .. code-block:: console
  
            python3 -m venv sphinx-env
            source sphinx-env/bin/activate

      .. group-tab:: Windows

         .. code-block:: console
    
            python3 -m venv sphinx-env
            .\sphinx-env\Scripts\activate

   Then confirm you're in the virtual environment:
   
   .. tabs::

      .. group-tab:: Unix/macOS
      
         .. code-block:: console
  
            which python

      .. group-tab:: Windows

         .. code-block:: console
    
            where python

   To leave the virtual environment:
   
    .. code-block:: console

       deactivate

3. Install `Sphinx <https://www.sphinx-doc.org>`_, `sphinxcontrib-matlabdomain <https://pypi.org/project/sphinxcontrib-matlabdomain/>`_, `sphinx-tabs <https://pypi.org/project/sphinx-tabs/>`_ and `sphinx-rtd-theme <https://pypi.org/project/sphinx-rtd-theme/>`_

   .. code-block:: console

      pip install -U sphinx
      pip install -U sphinxcontrib-matlabdomain
      pip install -U sphinx-tabs
      pip install -U sphinx-rtd-theme

4. Install `TeXLive <https://tug.org/texlive>`_ (for building LaTeX/PDF output).

5. In the :file:`docs/sphinx` directory, type:

    .. code-block:: console
    
       make html
       make latexpdf

   ... twice. That's right, you re-run ``make html`` after the LaTeX build so it can pick up the links to the PDF and then re-run ``make latexpdf`` to ensure that tables of contents, cross-references, etc. are all up to date.
   
   If everything builds properly, you should find the PDF manuals in :file:`docs/sphinx/build/latex` and the HTML documentation in :file:`docs/sphinx/build/html` (start with :file:`index.html`).
