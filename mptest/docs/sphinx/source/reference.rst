.. _sec_mptest_reference:

Reference
=========


Testing Functions
-----------------

Use these functions in implementing your own tests.

.. toctree::

   functions/t_begin
   functions/t_end
   functions/t_file_match
   functions/t_is
   functions/t_ok
   functions/t_run_tests
   functions/t_skip
   functions/t_str_match


Other Functions
---------------

Use these functions to test for availability and version information for
optional functionality and to check the version of the installed MP-Test.

.. toctree::

   functions/have_feature
   functions/mptestver


Tests of MP-Test
----------------

These functions test that MP-Test is installed and functioning as expected.

.. toctree::

   functions/test_mptest
   functions/t_have_feature
   functions/t_test_fcns


Private Functions
-----------------

The following are private functions that implement detection of specific
optional functionality. They are not intended to be called directly, but
rather are used to extend the capabilities of :func:`have_feature` *(see above)*.

.. toctree::

   functions/have_feature_matlab
   functions/have_feature_octave
