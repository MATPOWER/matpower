Introduction
============

The purpose of this *Developer's Manual* is to provide an understanding of the internal design of |MATPOWER| for users who wish to help with the development of |MATPOWER| or for those who would like to customize, modify or add to the functionality of |MATPOWER| in any way.

The |MATPOWER-Users-Manual|, on the other hand, is your starting point if you simply want to use |MATPOWER| without modification or customization.

For reference documentation on each class and function in |MATPOWER|, see the
|MATPOWER-Ref-Manual|.


Development Environment
-----------------------

|MATPOWER| is implemented in the Matlab language, designed for scientific computing. It requires either |MATLAB(R)|, a commercial product from `The MathWorks <https://www.mathworks.com/>`_, or the free, open-source `GNU Octave <https://www.octave.org>`_ to run.

|MATPOWER| and its related software packages are developed as open-source projects on GitHub under the `MATPOWER Development <https://github.com/MATPOWER>`_ GitHub organization. Some projects are included in others using `git subrepo <https://github.com/ingydotnet/git-subrepo>`_.

:numref:`tab_github_repos` provides an overview of the various repositories and their relationships to each other. Note that the main |gh-matpower|_ repository contains all of the others as subrepos, except for |gh-matpower-extras|_, which is, however, included when you download the ZIP file for a numbered |MATPOWER| release. 

.. _tab_github_repos:
.. list-table:: |MATPOWER| GitHub Repositories
   :widths: 22 78
   :header-rows: 1
   :class: longtable

   * - Repository
     - Description
   * - |gh-matpower|_
     - Main |MATPOWER| repository. Depends on |gh-mptest|_ , |gh-mips|_, and |gh-mp-opt-model|_, which are included as subrepos, along with |gh-most|_, and |gh-mp-docs-shared|_.
   * - |gh-mptest|_
     - Functions for implementing unit testing in MATLAB or Octave, with generalized mechanism for testing for optional functionality and corresponding versions, i.e. :func:`have_feature`. Required by all of the other projects.
   * - |gh-mips|_
     - |MIPSname| (|MIPS|), a nonlinear primal-dual interior point solver used as the default solver for AC OPF problems. Also includes a wrapper function for several linear equation solvers. Depends on |gh-mptest|_.
   * - |gh-mp-opt-model|_
     - |MPOM|, an easy-to-use, object-oriented interface for building and solving mathematical programming and optimization problems. Also includes a unified interface for calling numerous LP, QP, mixed-integer and nonlinear solvers, with the ability to switch solvers simply by changing an input option. Depends on |gh-mptest|_ and |gh-mips|_.
   * - |gh-most|_
     - |MOSTname| (|MOST|), a framework for solving generalized steady-state electric power scheduling problems. Depends on |gh-mptest|_, |gh-mp-opt-model|_ and |gh-matpower|_.
   * - |gh-matpower-extras|_
     - |MATPOWER| Extras, a collection of contributed and/or unsupported |MATPOWER|-related functions and packages. Note that some of the extras have their own separate repositories and are actually included here as subrepos. Depends on |gh-mptest|_ and |gh-matpower|_.
   * - |gh-mp-docs-shared|_
     - Defines common resources used for the Sphinx documentation and included as a subrepo in :file:`docs/sphinx/source` in all of the projects.

In general, each repository has two permanent branches, ``master`` and ``release``, where ``release`` points to the latest stable release and ``master`` contains any unreleased but hopefully stable updates. Each numbered release also has an associated git tag.


Conventions
-----------

Because |MATPOWER| is intended to run unmodified on either |MATLAB| or GNU Octave, it is important to stick to syntax and functionality that are supported by both.

We use :ml:`classdef` syntax supported by both to define classes and, in methods, we use ``obj`` as the variable name representing the object. Most of the classes are defined in the ``mp`` package/namespace.

All classes, methods, properties, and functions include a help section that can be accessed by the ``help`` and ``doc`` commands *and* processed by Sphinx to produce HTML and PDF reference documentation. For a class, it summarizes the purpose and overall functionality provided by the class along with lists of the properties and methods. For a function or method, it describes the inputs, outputs and what the function or method does. The :func:`run_mp` function and the :class:`mp.task` class provide examples of this reference documentation. *Hint: Click the GitHub icon in the upper right corner of the reference manual page to see the source.*

All functionality should be covered by at least one of the automated tests.

See the |Contrib_Guide|_ for more information on contributing to the |MATPOWER| project.


.. |gh-matpower| replace:: **matpower**
.. _gh-matpower: https://github.com/MATPOWER/matpower
.. |gh-mp-opt-model| replace:: **mp-opt-model**
.. _gh-mp-opt-model: https://github.com/MATPOWER/mp-opt-model
.. |gh-mips| replace:: **mips**
.. _gh-mips: https://github.com/MATPOWER/mips
.. |gh-most| replace:: **most**
.. _gh-most: https://github.com/MATPOWER/most
.. |gh-mptest| replace:: **mptest**
.. _gh-mptest: https://github.com/MATPOWER/mptest
.. |gh-matpower-extras| replace:: **matpower-extras**
.. _gh-matpower-extras: https://github.com/MATPOWER/matpower-extras
.. |gh-mp-docs-shared| replace:: **mp-docs-shared**
.. _gh-mp-docs-shared: https://github.com/MATPOWER/mp-docs-shared
.. |Contrib_Guide| replace:: |MATPOWER| Contributors Guide
.. _Contrib_Guide: https://github.com/MATPOWER/matpower/blob/master/CONTRIBUTING.md
