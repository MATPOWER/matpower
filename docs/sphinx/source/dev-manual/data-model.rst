.. _sec_data_model:

Data Model Object
=================

.. toctree::
   :maxdepth: 2

The data model is essentially the internal representation of the input data provided by the user for the given simulation or optimization run and the output presented back to the user upon completion. It corresponds roughly to the :ml:`mpc` (|MATPOWER| case) and :ml:`results` structs used throughout the legacy |MATPOWER| implementation, but encapsulated in an object with additional functionality. It includes tables of data for each type of element in the system.


Data Models
-----------

A data model object is primarily a container for data model element objects. All data model classes inherit from :class:`mp.data_model` and therefore also from :class:`mp.element_container`, and may be task-specific, as shown in :numref:`fig_data_model_classes`. For a simple power flow problem, :class:`mp.data_model` is used directly as a concrete class. For CPF and OPF problems, subclasses are used. In the case of CPF, :class:`mp.data_model_cpf` encapsulates both the base and the target cases. In the case of the OPF, :class:`mp.data_model_opf` includes additional input data, such as generator costs, and output data, such as nodal prices and shadow prices on line flow contraints.

.. _fig_data_model_classes:
.. figure:: figures/data-model-classes.*
   :alt: Data Model Classes
   :align: center
   :width: 300px

   Data Model Classes

By convention, data model variables are named :ml:`dm` and data model classes begin with :ml:`mp.data_model`.


Building a Data Model
^^^^^^^^^^^^^^^^^^^^^

There are two steps to building a data model. The first is to call the constructor of the desired data model class, without arguments. This initializes the :attr:`element_classes` property with a list of data model element classes. This list can be modified before the second step, which is to call the :meth:`build` method, passing in the data and a corresponding data model converter object.

.. _code_data_model_build:
.. code-block::

   dmc = mp.dm_converter_mpc2().build();
   dm = mp.data_model();
   dm.build('case9', dmc);

The :meth:`build` method proceeds through the following stages sequentially, looping through each element at every stage.

   1. **Create** – Instantiate each element object and add it to the :attr:`elements` property of the :ml:`dm`.
   2. **Import** – Use the corresponding data model converter element to read the data into each element's table(s).
   3. **Count** – Determine the number of instances of each element present in the data, store it in the element's :attr:`nr` property, and remove the element type from :attr:`elements` if the count is 0.
   4. **Initialize** – Initialize the (online/offline) status of each element and create a mapping of ID to row index in the :attr:`ID2i` element property.
   5. **Update status** – Update status of each element based on connectivity or other criteria and define element properties containing number and row indices of online elements (:attr:`n` and :attr:`on`), indices of offline elements (:attr:`off`), and mapping (:attr:`i2on`) of row indices to corresponding entries in :attr:`on` or :attr:`off`.
   6. **Build parameters** – Extract/convert/calculate parameters as necessary for online elements from the original data tables (e.g. p.u. conversion, initial state, etc.) and store them in element-specific properties.


System Level Parameters
^^^^^^^^^^^^^^^^^^^^^^^

There are a few system level parameters such as the system per-unit power base that are stored in data model properties. Balanced single-phase model elements, typical in transmission systems, use an MVA base found in :attr:`base_mva`. Unbalanced three-phase model elements, typical in distribution systems, use a kVA base found in :attr:`base_kva`. Models with both types of elements, therefore, use both properties.


Printing a Data Model
^^^^^^^^^^^^^^^^^^^^^

The :class:`mp.data_model` provides a :meth:`pretty_print` method for displaying the model parameters to a pretty-printed text format. The result can be output either to the console or to a file.

The output is organized into sections and each element type controls its own output for each section. The default sections are:

- **cnt** - count, number of online, offline, and total elements of this type
- **sum** - summary, e.g. total amount of capacity, load, line loss, etc.
- **ext** - extremes, e.g. min and max voltages, nodal prices, etc.
- **det** - details, table of detailed data, e.g. voltages, prices for buses, dispatch, limits for generators, etc.


Data Model Elements
-------------------

A data model element object encapsulates all of the input and output data for a particular element type. All data model element classes inherit from :class:`mp.dm_element` and each element type typically implements its own subclass. A given data model element object contains the data for all instances of that element type, stored in one or more *table* data structures. [#]_ So, for example, the data model element for generators contains a table with the generator data for all generators in the system, where each table row corresponds to an individual generator.

By convention, data model element variables are named :ml:`dme` and data model element classes begin with :ml:`dme`. :numref:`fig_data_model_element_classes` shows the inheritance relationships between a few example data model element classes. Here the :class:`mp.dme_bus`, :class:`mp.dme_gen` and :class:`mp.dme_load` classes are used for PF and CPF runs, while the OPF requires task-specific subclasses of each.

.. _fig_data_model_element_classes:
.. figure:: figures/data-model-element-classes.*
   :alt: Data Model Element Classes
   :align: center
   :width: 450px

   Data Model Element Classes

Data Tables
^^^^^^^^^^^

Typically, a data model element has at least one main table, stored in the :attr:`tab` property. Each row in the table corresponds to an individual element and the columns are the parameters. In general, |MATPOWER| attempts to follow the parameter naming conventions outlined in *The Common Electric Power Transmission System Model* (CTM) [CTM]_. The following parameters (table columns) are shared by all data model elements.

  - **uid** – positive integer serving as a unique identifier for the element
  - **name** – optional string identifier for the element
  - **status** – 0 or 1, on/off-line status of the element
  - **source_uid** – implementation specific *(e.g. sometimes used to map back to a specific record in the source data)*

Properties
^^^^^^^^^^

The table below includes additional properties, besides the main table :attr:`tab`, found in all data model elements.

.. list-table:: 
   :widths: 12 88
   :header-rows: 1
   :class: longtable

   * - Property
     - Description
   * - :attr:`nr`
     - number of rows in the table, i.e. the *total* number of elements of the type
   * - :attr:`n`
     - number of *online* elements of the type
   * - :attr:`on`
     - vector of row indices of online elements
   * - :attr:`off`
     - vector of row indices of offline elements
   * - :attr:`ID2i`
     - :math:`M \times 1` vector mapping IDs to row indices, where :math:`M` is the largest ID value
   * - :attr:`i2on`
     - :math:`n_r \times 1` vector mapping row indices to the corresponding index into the :attr:`on`
       vector *(for online elements)* or :attr:`off` vector *(for offline elements)*
   * - :attr:`tab`
     - main data table

Methods
^^^^^^^

A data model element also has a :meth:`name` method that returns the name of the element type under which it is entered in the data model (container) object. For example, the name returned for the :class:`mp.dme_gen` class is :ml:`'gen'`, which means the corresponding data model element object is found in :ml:`dm.elements.gen`.

There are also methods named :meth:`label` and :meth:`labels` which return user visible names for singular and plural instances of the element used when pretty-printing. For :class:`mp.dme_gen`, for example, these return :ml:`'Generator'` and :ml:`'Generators'`, respectively.

.. note::

   Given that these name and label methods simply return a character array, it might seem logical to implement them as Constant properties instead of methods, but this would prevent the value from being overridden by a subclass, in effect precluding the option to create a new element type that inherits from an existing one.

The :meth:`main_table_var_names` method returns a cell array of variable names defining the columns of the main data table. These are used by the corresponding data model converter element to import the data.

There are also methods that correspond to the build steps for the data model container object, :meth:`count`, :meth:`initialize`, :meth:`init_status`, :meth:`update_status`, and :meth:`build_params`, as well as those for pretty printing output, :meth:`pretty_print`, etc.

.. [#] Implemented using the built-in :class:`table` and included :class:`mp_table` classes, respectively, under |MATLAB| and GNU Octave. See also :func:`mp_table_class`.

.. [CTM] Carleton Coffrin, et. al., "The Common Electric Power Transmission System Model," *work in progress*. Available at: https://www.overleaf.com/project/5d94e3765cb3ba000129df3c.

