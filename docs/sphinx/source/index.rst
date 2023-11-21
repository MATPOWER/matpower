.. MP-Docs documentation master file, created by
   sphinx-quickstart on Mon Jul 11 17:08:23 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

    .. image:: MATPOWER-md.png

########################
|MATPOWER| Documentation
########################

..
    .. note::

       This is a very incomplete **work-in-progress**!

.. toctree::
   :maxdepth: 1

   Home (matpower.org) <https://matpower.org>

------------------------------------

Overview
========

The documentation for |MATPOWER| and its related packages are distributed across several manuals and How To Guides. The goal is to make all of this documentation available in HTML format, with some of the manuals available in PDF format as well. To facilitate this goal, the source for new documentation will be written in `reStructured Text <https://www.writethedocs.org/guide/writing/reStructuredText/>`_ format for processing by `Sphinx <https://www.sphinx-doc.org/>`_ to build the HTML and PDF versions from a single source.

Currently, most of the User's Manuals are still available only as PDF, built from the legacy LaTeX source.

|MATPOWER| is a package of free, open-source Matlab-language M-files for solving steady-state power system simulation and optimization problems. It is built on |MPOM>|, a package for constructing and solving mathematical programming and optimization problems. |MPOM>| in turn relies on |MIPS>| (|MIPSname|) as a default solver. And |MOST>| (|MOSTname|), built on top of |MATPOWER|, is a framework for solving generalized steady-state electric power scheduling problems. All of them rely on |MPTEST>| as a software testing framework.

------------------------------------

.. toctree::
   :maxdepth: 1
   :caption: Manuals

   users-manual/index
   dev-manual/index
   ref-manual/index

.. toctree::
   :maxdepth: 1

   mptest/index
   mips/index
   mpom/index
   most/index

------------------------------------

.. toctree::
   :maxdepth: 2
   :caption: How To Guides

   howto/element
   howto/constraint
   howto/extension
   howto/builddocs

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

       users-manual-legacy/index
       _installation
       _reference
       _api


    Indices and tables
    ==================

    * :ref:`genindex`
    * :ref:`modindex`
    * :ref:`search`
