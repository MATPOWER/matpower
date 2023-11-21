########################
|MATPOWER| User's Manual
########################

.. only:: html

   .. image:: ../MATPOWER-md.png

.. note::

   The new web-based version of the User's Manual is not yet available. Please, continue to use the PDF version of the |MUM| for now.


Legacy Framework
================

.. note::

   User documentation for the **legacy framework** is found in the legacy |MUM|.

The legacy |MUM| decribes the **legacy** |MATPOWER| framework and interface, that is, the functionality and features fully compatible with |MATPOWER| 7.x or earlier.

For example, you can run an AC OPF on the nine bus case as follows.

.. code-block::

   runopf('case9')

By default, this uses the **legacy** |MATPOWER| framework, for full compatibility with legacy features and extensions, with MP-Core providing the underlying modeling. However, MP-Core can be bypassed completely, with legacy code being used for everything by passing an option to a specific run.

.. code-block::

   mpopt = mpoption('exp.use_legacy_core', 1);
   runopf('case9', mpopt)

Or MP-Core can be disabled globally for the current session with the following command, in which case no legacy commands will use MP-Core internally.

.. code-block::

   have_feature('mp_core', 0)


New Flexible Framework
======================

.. note::

   User documentation for the new **flexible framework** is found mainly in the |MATPOWER-Ref-Manual|, with some additional information in the |MATPOWER-Dev-Manual|.

To run an AC OPF on the nine bus case in this mode you need to use the version of the ``run`` commands with an underscore in the name.

.. code-block::

   run_opf('case9')

See the :ref:`ref_top_level_functions` section in the |MATPOWER-Ref-Manual| to get started.
