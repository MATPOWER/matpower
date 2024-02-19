.. _sec_customizing:

Customizing |MATPOWER|
======================

With the *object-oriented* |/MATPOWER/| *core architecture* and its explicit three layer modeling outlined in :numref:`sec_architecture`, the flexibility and customizability of |MATPOWER| has increased dramatically. 

New functionality can be added or existing functionality modified simply by adding new classes and/or subclassing existing ones. This approach can be used to add or modify elements, problem formulations, and tasks.


Default Class Determination
---------------------------

In order to customize the behavior it is important to understand how |MATPOWER| selects which classes to instantiate when running a particular task. There are default specifications for each of the various types of objects, as well as several ways to override those defaults. The default, described below, is illustrated in :numref:`fig_class_selection`.

.. _fig_class_selection:
.. figure:: figures/class-selection.*
   :alt: Determination of Default Classes
   :align: center

   Determination of Default Classes


Task Class
^^^^^^^^^^

First of all, at the top level, the **task class** is specified directly by the user through the function used to invoke the run. In fact, :func:`run_pf`, :func:`run_cpf`, and :func:`run_opf` are simple one-line wrappers around the :func:`run_mp` function. The only difference between the three is the value of the ``task_class`` argument, a handle to the corresponding task constructor, passed into :func:`run_mp`.

This means that a new task class can be used simply by invoking :func:`run_mp`, either directly or via a new wrapper, with the task constructor as the ``task_class`` argument.


Model and Data Converter Classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The task class has methods that determine the classes used to create the data model converter and the three model objects. For each, there are two methods involved in determining the specific class to use, a *main* method and a *default* method. The *main* method calls the *default* method to get a handle to the constructor for the task's default class, but then allows that value to be overridden by |MATPOWER| extensions or |MATPOWER| options.

.. _tab_task_methods:
.. list-table:: Class Specification Methods of a Task
   :widths: 40 60
   :header-rows: 1
   :class: longtable

   * - Method
     - Description
   * - :meth:`dm_converter_class() <mp.task.dm_converter_class>`
     - Returns the final class for the data model converter, after any overrides of the default.
   * - :meth:`data_model_class() <mp.task.data_model_class>`
     - Returns the final class for the data model, after any overrides of the default.
   * - :meth:`network_model_class() <mp.task.network_model_class>`
     - Returns the final class for the network model, after any overrides of the default.
   * - :meth:`math_model_class() <mp.task.math_model_class>`
     - Returns the final class for the math model, after any overrides of the default.
   * - :meth:`dm_converter_class_mpc2_default() <mp.task.dm_converter_class_mpc2_default>`
     - Returns the *default* class for the data model converter for this task. Note that this is specific to the data format. Each data format would have it's own "default" method.
   * - :meth:`data_model_class_default() <mp.task.data_model_class_default>`
     - Returns the *default* class for the data model for this task.
   * - :meth:`network_model_class_default() <mp.task.network_model_class_default>`
     - Returns the *default* class for the network model for this task.
   * - :meth:`math_model_class_default() <mp.task.math_model_class_default>`
     - Returns the *default* class for the math model for this task.

:numref:`tab_task_methods` shows the methods that determine the classes for each of the 4 objects. Each method returns a handle to a class constructor. In general, the *main* methods (the first 4 in the table) are inherited from :class:`mp.task` and only the *default* methods (the last 4) would be overridden to customize a task with new model or data model converter classes.


.. _sec_element_classes:

Element Classes
^^^^^^^^^^^^^^^

Each of the element container objects, that is the data model converter and the 3 model objects, contains a set of *elements*. The classes used to construct these elements are determined by the container class. Each container class inherits from :class:`mp.element_container`, and as such it has an :attr:`element_classes <mp.element_container.element_classes>` property, which is a cell array populated by the  constructor with handles to constructors for the elements. This means that a container subclass can, by overriding its constructor, modify the list of element classes provided by its parent.

The elements are instantiated by a call to the container object's :meth:`build` method, so the resulting set can be customized at runtime by modifying the list in :attr:`element_classes <mp.element_container.element_classes>` after the container object is created and before its :meth:`build` method is called.

This is done using **element class modifiers**, specified either by |MATPOWER| extensions or |MATPOWER| options. There are 3 types of element class modifiers, for adding, deleting or replacing an entry in an :attr:`element_classes <mp.element_container.element_classes>` property. The 3 types are described in :numref:`tab_element_class_modifiers`.


.. _tab_element_class_modifiers:
.. list-table:: Element Class Modifiers
   :widths: 10 33 57
   :header-rows: 1
   :class: longtable

   * - Action
     - Value
     - Description
   * - **add**
     - ``@new_class``
     - Appends ``@new_class`` to the end of the list.
   * - **delete**
     - ``'old_class'``
     - For each element ``E`` in the list, if :ml:`isa(E(), 'old_class')` is true, element ``E`` is deleted from the list.
   * - **replace**
     - ``{@new_class, 'old_class'}``
     - For each element ``E`` in the list, if :ml:`isa(E(), 'old_class')` is true, element ``E`` is replaced with ``@new_class``.

Typically, multiple element class modifiers can be provided in a cell array and they are processed sequentially to modify the existing list by the :meth:`modify_element_classes() <mp.element_container.modify_element_classes>` from :class:`mp.element_container`.


Customization via |MATPOWER| Options
------------------------------------

In addition to the |MATPOWER| options previously available that affect the formulation of the problem (e.g. polar vs. cartesian voltage representation, or current vs. power balance), there are several experimental options that can be used to directly modify the classes coming from the default class selection process outlined above. These options, summarized in :numref:`tab_custom_class_options`, are specified by assigning them directly to an existing |MATPOWER| options struct ``mpopt`` as optional fields in ``mpopt.exp``. They must be assigned directly, since :func:`mpoption` does not recognize them.

.. _tab_custom_class_options:
.. list-table:: Class Customization Options
   :widths: 25 75
   :header-rows: 1
   :class: longtable

   * - Option
     - Description
   * - ``dm_converter_class``
     - function handle for data model converter constructor
   * - ``data_model_class``
     - function handle for data model constructor
   * - ``network_model_class``
     - function handle for network model constructor
   * - ``math_model_class``
     - function handle for math model constructor
   * - ``dmc_element_classes``
     - element class modifier(s) [#]_ for data model converter elements
   * - ``dm_element_classes``
     - element class modifier(s) [1]_ for data model elements
   * - ``nm_element_classes``
     - element class modifier(s) [1]_ for network model elements
   * - ``mm_element_classes``
     - element class modifier(s) [1]_ for math model elements
   * - ``exclude_elements``
     - cell array of names of elements to exclude from all four container objects, i.e. char arrays that match the :attr:`name` property of the element(s) to be excluded


.. _sec_extensions:

|MATPOWER| Extensions
---------------------

The *flexible* |/MATPOWER/| *framework* summarized in :numref:`sec_two_frameworks` introduces a |*MATPOWER*| **extension** API as a way to bundle a set of class additions and modifications together into a single named package.

For example, the :class:`mp.xt_reserves` class and those it references, adds co-optimization of fixed zonal reserves to the standard OPF problem, as previously implemented by :func:`toggle_reserves` and :func:`run_opf_w_res` in |MATPOWER| 7.1 and earlier using its legacy OPF callback functions. To invoke an OPF with the :class:`mp.xt_reserves` extension, you simply pass the extension object as an optional argument into the :func:`run_opf` function.

.. code-block::

   run_opf(mpc, mpopt, 'mpx', mp.xt_reserves);

A |MATPOWER| extension is a subclass of :class:`mp.extension`, which implements a very simple interface consisting of nine methods. Five of them return a single class constructor handle, and the other four return a cell array of element class modifiers, described above in :numref:`tab_element_class_modifiers`.

The methods are summarized in :numref:`tab_ext_methods`

.. _tab_ext_methods:
.. list-table:: |MATPOWER| Extension Methods
   :widths: 25 75
   :header-rows: 1
   :class: longtable

   * - Method
     - Description
   * - :meth:`task_class() <mp.extension.task_class>`
     - Returns a handle to the constructor for the task object.
   * - :meth:`dm_converter_class() <mp.extension.dm_converter_class>`
     - Returns a handle to the constructor for the data model converter.
   * - :meth:`data_model_class() <mp.extension.data_model_class>`
     - Returns a handle to the constructor for the data model.
   * - :meth:`network_model_class() <mp.extension.network_model_class>`
     - Returns a handle to the constructor for the network model.
   * - :meth:`math_model_class() <mp.extension.math_model_class>`
     - Returns a handle to the constructor for the math model.
   * - :meth:`dmc_element_classes() <mp.extension.dmc_element_classes>`
     - Returns a cell array of element class modifiers for data model converter elements.
   * - :meth:`dm_element_classes() <mp.extension.dm_element_classes>`
     - Returns a cell array of element class modifiers for data model elements.
   * - :meth:`nm_element_classes() <mp.extension.nm_element_classes>`
     - Returns a cell array of element class modifiers for network model elements.
   * - :meth:`mm_element_classes() <mp.extension.mm_element_classes>`
     - Returns a cell array of element class modifiers for math model elements.

Even something as complex as adding three-phase unbalanced buses, lines, loads and generators for multiple formulations of PF, CPF, and OPF problems can be implemented in terms of a single |MATPOWER| extension. Please see :class:`mp.xt_3p` for an example.


.. [#] Either a single element class modifier or a cell array of element class modifiers.

..
    Careful the footnote above is explicitly numbered as [1]_ in several
    references above (to avoid repeating the footnote itself).
