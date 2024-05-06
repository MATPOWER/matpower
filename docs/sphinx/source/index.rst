.. MP-Docs documentation master file, created by
   sphinx-quickstart on Mon Jul 11 17:08:23 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

    .. image:: MATPOWER-md.png

########################
|MATPOWER| Documentation
########################

.. toctree::
   :caption: Home

   matpower.org <https://matpower.org>

.. toctree::
   :caption: Get Started

   Get Started <https://matpower.org/about/get-started/>

------------------------------------

Overview
========

..
    |MATPOWER| is a package of free, open-source Matlab-language M-files for solving steady-state power system simulation and optimization problems. It is built on |MPOM>|, a package for constructing and solving mathematical programming and optimization problems. |MPOM>| in turn relies on |MIPS>| (|MIPSname|) as a default solver. And |MOST>| (|MOSTname|), built on top of |MATPOWER|, is a framework for solving generalized steady-state electric power scheduling problems. All of them rely on |MPTEST>| as a software testing framework.


|MATPOWER| is a package of free, open-source Matlab-language functions and classes for simulating and optimizing steady-state power system problems. It includes the following sub-packages.

.. list-table::
   :widths: 22 58 20
   :header-rows: 1
   :class: longtable

   * - Package
     - Description
     - Requires
   * - |MATPOWER>|
     - steady-state power system simulation and optimization
     - |MPOM|
   * - |MPOM>|
     - tools for constructing and solving mathematical programming and optimization problems
     - |MIPS|
   * - |MIPS>|
     - |MIPSname|, nonlinear programming (NLP) solver
     - |MPTEST|
   * - |MPTEST>|
     - software testing framework and utility functions
     - |MATLAB|/Octave
   * - |MOST>|
     - |MOSTname|, framework for generalized steady-state electric power scheduling problems
     - |MATPOWER|
   * - |MATPOWER-Extras>|
     - additional contributed packages for use with |MATPOWER|
     - |MATPOWER|

The documentation for |MATPOWER| can be found in the various Manuals, How To Guides, and Tech Notes on this site. Each function or class has its own reference documentation which can be accessed via the built-in ``help`` command in |MATLAB| or Octave, or by consulting the corresponding Reference manual.

Changes to |MATPOWER| in each released version are summarized in the `release notes <https://github.com/MATPOWER/matpower/blob/master/docs/relnotes>`_, found on GitHub in ``docs/relnotes`` and in Appendix H of the |MUM|. A complete, detailed change log, even for unreleased versions, is available in the `CHANGES.md <https://github.com/MATPOWER/matpower/blob/master/CHANGES.md>`_ file.


.. note::

   The goal is to make all |MATPOWER| documentation available in HTML format, with some of the manuals also available as PDF. To facilitate this, new documentation is written in `reStructured Text <https://www.writethedocs.org/guide/writing/reStructuredText/>`_ format for processing by `Sphinx <https://www.sphinx-doc.org/>`_, so the HTML and PDF versions can be built from a single source.

   Currently, most of the User's Manuals are still only available as PDF, built from the legacy LaTeX source.


------------------------------------

.. toctree::
   :maxdepth: 1
   :caption: Manuals

   users-manual/index
   dev-manual/index
   ref-manual/index

More Manuals
------------

.. toctree::
   :maxdepth: 1

   mptest/index
   mips/index
   mpom/index
   most/index

All Legacy PDF Manuals
----------------------

.. toctree::

   All Legacy PDF Manuals <https://matpower.org/doc/manuals/>

------------------------------------

.. toctree::
   :maxdepth: 1
   :caption: How To Guides

   howto/element
   howto/add-constraint
   howto/extension
   howto/three-phase
   howto/builddocs

------------------------------------

.. toctree::
   :maxdepth: 2
   :caption: Tech Notes

   tech-notes

------------------------------------

Publications
============

For additional |MATPOWER|-related publications, see:

- `Publications <https://matpower.org/publications>`_

Please `cite <https://matpower.org/citing>`_ |MATPOWER| in your own publications derived from the use of |MATPOWER|.

- `Citing <https://matpower.org/citing>`_ |MATPOWER|

------------------------------------

.. toctree::
   :caption: Other Links

   Donate <https://matpower.org/sponsor>
   Downloads <https://matpower.org/download>
   GitHub Project <http://github.com/MATPOWER/matpower>


..
    ------------------------------------

    And other junk

    .. toctree::
       :maxdepth: 2
       :caption: Other Junk

       howto/constraint
       users-manual-legacy/index
       _installation
       _reference
       _api


    Indices and tables
    ==================

    * :ref:`genindex`
    * :ref:`modindex`
    * :ref:`search`
