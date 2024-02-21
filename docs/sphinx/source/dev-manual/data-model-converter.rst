.. _sec_dm_converter:

Data Model Converter Object
===========================

A data model converter provides the ability to convert data between a data model and a specific data source or format, such as the PSS/E RAW format or version 2 of the |MATPOWER| case format. It is used, for example, during the **import** stage of the data model build process.

Data Model Converters
---------------------

A data model converter object is primarily a container for data model converter element objects. All data model converter classes inherit from :class:`mp.dm_converter` and therefore also from :class:`mp.element_container` and they are specific to the type or format of the data source, as shown in :numref:`fig_dm_converter_classes`. In this example, the PSS/E RAW format converter has not yet been implemented, but is shown here for illustration.

.. _fig_dm_converter_classes:
.. figure:: figures/dm-converter-classes.*
   :alt: Data Model Converter Classes
   :align: center
   :width: 400px

   Data Model Converter Classes

By convention, data model converter variables are named ``dmc`` and data model converter class names begin with ``mp.dm_converter``.

Building a Data Model Converter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A data model converter object is created in two steps. The first is to call the constructor of the desired data model converter class, without arguments. This initializes the :attr:`element_classes <mp.element_container.element_classes>` property with a list of data model converter element classes. This list can be modified before the second step, which is to call the :meth:`build() <mp.dm_converter.build>` method, also without parameters, which simply instantiates and adds the set of element objects indicated in :attr:`element_classes <mp.element_container.element_classes>`. Once it has been created, it is ready to be used for its two primary functions, namely *import* and *export*.

.. _code_data_model_build:
.. code-block::

   dmc = mp.dm_converter_mpc2();
   dmc.build();

Importing Data
^^^^^^^^^^^^^^

The :meth:`import() <mp.dm_converter.import>` method is called automatically by the :meth:`build() <mp.data_model.build>` method of the data model object. It takes a data model object and a data source and updates the data model by looping through its element objects and calling each element's own :meth:`import() <mp.dmc_element.import>` method to import the element's data from the data source into the corresponding data model element. For a |MATPOWER| case struct it would like like this.

.. _code_dmc_import:
.. code-block::

   mpc = loadcase('case9');
   dm = dmc.import(dm, mpc);

Exporting Data
^^^^^^^^^^^^^^

Conversely, the :meth:`export() <mp.dm_converter.export>` method takes the same inputs but returns an updated data source, once again looping through its element objects and calling each element's own :meth:`export() <mp.dmc_element.export>` method to export data from the corresponding data model element to the respective portion of the data source.

.. _code_dmc_export:
.. code-block::

   mpc = dmc.export(dm, mpc);

Calling :meth:`export() <mp.dm_converter.export>` without passing in a data source will initialize one from scratch.

.. _code_dmc_export_init:
.. code-block::

   mpc = dmc.export(dm);


.. _sec_dmc_element:

Data Model Converter Elements
-----------------------------

A data model converter element object implements the functionality needed to import and export a particular element type from and to a given data format. All data model converter element classes inherit from :class:`mp.dmc_element` and each element type typically implements its own subclass.

By convention, data model converter element variables are named ``dmce`` and data model converter element class names begin with ``mp.dmce``. :numref:`fig_dm_converter_classes` shows the inheritance relationships between a few example data model converter element classes. Here the PSS/E classes have not yet been implemented, but are shown here for illustration.

.. _fig_dm_converter_element_classes:
.. figure:: figures/dm-converter-element-classes.*
   :alt: Data Model Converter Element Classes
   :align: center
   :width: 600px

   Data Model Converter Element Classes

Data Import Specifications
^^^^^^^^^^^^^^^^^^^^^^^^^^

The default :meth:`import() <mp.dmc_element.import>` method for a data model converter element first calls the :meth:`get_import_spec() <mp.dmc_element.get_import_spec>` method to get a struct containing the specifications that define the details of the import process. This specification is then passed to :meth:`import_table_values() <mp.dmc_element.import_table_values>` to import the data.

The import specifications include things like where to find the data in the data source, the number of rows, number of columns, and possibly a row index vector for rows of interest, [#]_ and a map defining how to import each column of the main data table.

This map ``vmap`` is a struct returned by the :meth:`table_var_map() <mp.dmc_element.table_var_map>` method with fields matching the table column names for the corresponding data model element ``dme``. For example, if ``vn`` contains a variable, that is column, name, then :samp:`vmap.(vn) = {<value>}` defines how that data table column will be imported or initialized, as summarized in :numref:`tab_var_map` for different types of values.

.. _tab_var_map:
.. list-table:: Variable Map Values
   :widths: 25 75
   :header-rows: 1
   :class: longtable

   * - :samp:`{<value>}`
     - Description
   * - ``{'IDs'}``
     - Assign consecutive IDs starting at 1.
   * - :samp:`\\{'col', {c}\\}` *or*
   
       :samp:`\\{'col', {c}, {sf}\\}` *or*
   
       :samp:`\\{'col', {c}, {sf_fcn}\\}`
     - Copy the data directly from column :samp:`{c}` of data source, optionally scaling it by a numerical scale factor :samp:`{sf}`, or by the value returned by the function handle :samp:`{sf_fcn}`, called as :samp:`{sf_fcn(dmce, vn)}`.
   * - :samp:`\\{'cell', {val}\\}`
     - Create a cell array with each element initialized with :samp:`{val}`.
   * - :samp:`\\{'num', {n}\\}`
     - Create a numeric vector with each element initialized with numeric scalar :samp:`{n}`.
   * - :samp:`\\{'fcn', {ifn}\\}` *or*

       :samp:`\\{'fcn', {ifn}, {efn}\\}`
     - Assign the values returned by the import function handle in :samp:`{ifn}`, where the optional :samp:`{efn}` is the corresponding export function. The import and export functions are called as :samp:`{ifn(dmce, d, spec, vn)}` and :samp:`{efn(dmce, dme, d, spec, vn, ridx)}`, respectively, where :samp:`{d}` is the data source, :samp:`{spec}` is the import/export specification, and :samp:`{ridx}` is an optional vector of row indices.

The :meth:`table_var_map() <mp.dmc_element.table_var_map>` in :class:`mp.dmc_element` initializes each entry to ``{'col', []}`` by default, so subclasses only need to assign ``vmap.(vn){2}`` for columns that map directly from a column of the data source.


Data Export Specifications
^^^^^^^^^^^^^^^^^^^^^^^^^^

The default :meth:`export() <mp.dmc_element.export>` method first calls the :meth:`get_export_spec() <mp.dmc_element.get_export_spec>` method to get a struct containing the specifications that define the details of the export process. This specification is then passed to :meth:`export_table_values() <mp.dmc_element.export_table_values>` to export the data.

The export of data from a data model element back to the original data format is handled by the same variable map as the input, by default.

The :meth:`init_export_data() <mp.dmc_element.init_export_data>` method is used to initialize the relevant output data structure before exporting to it, if the :meth:`data_exists() <mp.dmc_element.data_exists>` method returns false.


.. [#] For example, when extracting loads from a bus matrix, where only certain buses have corresponding loads.