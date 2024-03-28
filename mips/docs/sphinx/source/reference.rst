.. _sec_mips_reference:

Reference
=========


Main MIPS Functions
-------------------

The MIPS nonlinear programming (NLP) solver, a linear solver API function, and a QP solver wrapper function for MIPS.

.. toctree::

   functions/mips
   functions/mplinsolve
   functions/qps_mips


Other Function
--------------

Use this function to check the version of the installed MIPS.

.. toctree::

   functions/mipsver


MIPS Examples
-------------

These are examples of using MIPS to solve a NLP.

.. toctree::

   functions/mips_example1
   functions/mips_example2


MIPS Tests
----------

These functions test that MIPS is installed and functioning as expected.

.. toctree::

   functions/test_mips
   functions/t_mips
   functions/t_mips_pardiso
   functions/t_mplinsolve
   functions/t_qps_mips


Private Functions
-----------------

The following are private functions that implement detection of specific
optional functionality. They are not intended to be called directly, but
rather are used to extend the capabilities of :func:`have_feature`.

.. toctree::

   functions/have_feature_lu_vec
   functions/have_feature_pardiso_legacy
   functions/have_feature_pardiso_object
   functions/have_feature_pardiso
