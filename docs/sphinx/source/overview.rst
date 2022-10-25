Overview
========

.. toctree::
   :maxdepth: 2

Overview of |MATPOWER| Documentation
------------------------------------

An overview of the |MATPOWER| Documentation site.


How To Compile the Documentation
--------------------------------

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

5. In the :file:`docs` directory, type:

    .. code-block:: console
    
       make html
       make latexpdf
