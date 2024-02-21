.. _sec_net_model:

Network Model Object
====================

The network model defines the states of and connections between network elements, as well as the parameters and functions defining the relationships between states and port injections. This network model with a unified structure is *the key* to a flexible modular design where model elements can simply define a few parameters and all of the mathematics involved in computing injections and their derivatives for a given state is handled automatically.

One of the unique features of the network model is that the network model object, which contains network model elements, **is a** network model element itself. Each network model element can optionally define the following:

- **nodes** to serve as network connection points
- **ports** for connecting to network nodes
- **states** which fully capture the element's operating condition

There are two types of states that make up the element's full state variable :math:`\X`, voltage states :math:`\V` associated with each port, and optional non-voltage states :math:`\Z`. The network model object inherits from :class:`mp_idx_manager` from |MPOM>| to track and index the nodes, ports, and states added by its elements, and the corresponding voltage and non-voltage state variables.

A given network model implements a specific network model **formulation**.
:numref:`fig_ac_net_element_model` shows the structure of an AC network model for an element with :math:`n_p` connection ports and :math:`n_z` non-voltage states.

.. _fig_ac_net_element_model:
.. figure:: figures/mp-element-ac-model.*
   :alt: AC Model for Element with n_p Ports
   :align: center
   :width: 360px

   AC Model for Element with :math:`n_p` Ports

The MP-Element technical note, |TN5| [TN5]_, includes a lot more detail on the network model and especially on the mathematics involved in the model formulations.


.. _sec_net_model_formulations:

Network Model Formulations
--------------------------

Each concrete network model element class, including the container class, inherits from a specific subclass of :class:`mp.form`. That is, it implements a specific network model formulation. For example, :numref:`fig_nm_form_classes` shows the corresponding classes for the three network model formulations currently implemented, (1) DC, (2) AC with cartesian voltage representation, and (3) AC with polar voltage representation.

.. _fig_nm_form_classes:
.. figure:: figures/nm-form-classes.*
   :alt: Network Model Formulation Classes
   :align: center
   :width: 380px

   Network Model Formulation Classes

By convention, network model formulation class names begin with ``mp.form``. It is the formulation class that defines the network model's parameters and methods for accessing them. It also defines the form of the state variables, real or complex, and methods for computing injections as a function of the state, and in the case of nonlinear formulations, corresponding derivatives as well.

All formulations share a common structure, illustrated in :numref:`fig_ac_net_element_model`, with ports, corresponding voltage states, non-voltage states, and functions of predefined form for computing port injections from the state.


.. _sec_nm_formulations_dc:

DC Formulation
^^^^^^^^^^^^^^

For the DC formulation, the state vector :math:`\x` is real valued and the port injection function is defined in terms of active power injections. The state begins with the :math:`n_p \times 1` vector :math:`\Va` of voltage angles at the :math:`n_p` ports, and may include an :math:`n_\z \times 1` real vector of additional state variables :math:`\z`, for a total of :math:`n_\x` state variables.

.. math::
   :label: eq_Xk_DC

   \x = \left[\begin{array}{c}
   \Va \\
   \z
   \end{array}\right]

The port injection function in this case defines the active power port injections as a linear function of a set of parameters :math:`\BB`, :math:`\KK` and :math:`\pv`, where :math:`\BB` is an :math:`n_p \times n_p` susceptance matrix, :math:`\KK` is an :math:`n_p \times n_\z` matrix coefficient for a linear power injection function, and :math:`\pv` is an :math:`n_p \times 1` constant power injection.

.. math::
   :label: eq_GkPx

   \gP(\x) &= \left[\begin{array}{cc}\BB & \KK\end{array}\right] \x + \pv \\
   &= \BB \Va + \KK \z + \pv

.. 
   \gP(\x) &= \left[\begin{array}{cc}\BB & \KK\end{array}\right] \x + \pv \\
   &= \BB \Va + \KK \z + \pv \\[2ex]
   \gP_\x & = \left[\begin{array}{cc}\BB & \KK\end{array}\right]


.. _sec_nm_formulations_ac:

AC Formulations
^^^^^^^^^^^^^^^

For the AC formulations, the state vector :math:`\X` is complex valued and there are two port injection functions, one for complex power injections and one for current injections, as shown in :numref:`fig_ac_net_element_model`. The state begins with the :math:`n_p \times 1` vector :math:`\V` of complex voltages at the :math:`n_p` ports, and may include an :math:`n_\Z \times 1` real vector of additional state variables :math:`\Z`, for a total of :math:`n_\X` state variables.

.. math::
   :label: eq_Xk_AC

   \X = \left[\begin{array}{c}
   \V \\
   \Z
   \end{array}\right]

The port injection functions for the model, both complex power injection :math:`\GS(\X)` and complex current injection :math:`\GI(\X)`, are defined by three terms, a linear current injection component :math:`\Ilin(\X)`, a linear power injection component :math:`\Slin(\X)`, and an arbitrary nonlinear component, :math:`\Snln(\X)` or :math:`\Inln(\X)`, respectively.

The linear current and power injection components are expressed in terms of the six parameters, :math:`\YY`, :math:`\LL`, :math:`\MM`, :math:`\NN`, :math:`\iv`, and :math:`\sv`. The admittance matrix :math:`\YY` and linear power coefficient matrix :math:`\MM` are :math:`n_p \times n_p`, linear coefficient matrices :math:`\LL` and :math:`\NN` are :math:`n_p \times n_\Z`, and :math:`\iv` and :math:`\sv` are :math:`n_p \times 1` vectors of constant current and power injections, respectively.

.. math::
   :label: eq_Ilin

   \Ilin(\X) &= \left[\begin{array}{cc}\YY & \LL\end{array}\right] \X + \iv \\
   &= \YY \V + \LL \Z + \iv

.. math::
   :label: eq_Slin

   \Slin(\X) &= \left[\begin{array}{cc}\MM & \NN\end{array}\right] \X + \sv \\
   &= \MM \V + \NN \Z + \sv

Note that the arbitrary *nonlinear* injection component, represented by either :math:`\Snln(\X)` or :math:`\Inln(\X)`, corresponds to a single set of injections represented either as a complex power injection or as a complex current injection, but not both. Since the functions represent the same set of injections, they are not additive components, but rather must be related to one another by the following relationship.

.. math::

   \Snln(\X) = \dV \conj{\left( \Inln(\X) \right)}

..
    We define :math:`\s(\X)` to be the power injection corresponding to the linear current term.

    .. math::
       :label: eq_SlinI

       \s(\X) = \dV \conj{\left( \Ilin(\X) \right)}

Complex Power Injections
''''''''''''''''''''''''

Then the port injection function for complex power can be written as follows.

.. math::
   :label: eq_GkS

   \GS(\X) &= \dV \conj{\left( \Ilin(\X) \right)} + \Slin(\X) + \Snln(\X) \\
   &= \dV \conj{\left( \YY \V + \LL \Z + \iv \right)} + \MM \V + \NN \Z + \sv + \Snln(\X)


Complex Current Injections
''''''''''''''''''''''''''

Similarly, the port injection function for complex current can be written as follows.

.. math::
   :label: eq_GkI

   \GI(\X) &= \Ilin(\X) + \cdiag{\Slin(\X)} \inVc + \Inln(\X) \\
   &= \YY \V + \LL \Z + \iv + \cdiag{\MM \V + \NN \Z + \sv} \inVc + \Inln(\X)

The derivatives of :math:`\Snln` and :math:`\Inln` are assumed to be provided explicitly, and the derivatives of the other terms of :eq:`eq_GkS` and :eq:`eq_GkI` are derived in [TN5]_.


Network Models
--------------

A network model object is primarily a container for network model element objects and *is itself* a network model element. All network model classes inherit from :class:`mp.net_model` and therefore also from :class:`mp.element_container`, :class:`mp_idx_manager`, and :class:`mp.nm_element`. Concrete network model classes are also formulation-specific, inheriting from a corresponding subclass of :class:`mp.form` as shown in :numref:`fig_net_model_classes`.

.. _fig_net_model_classes:
.. figure:: figures/net-model-classes.*
   :alt: Network Model Classes
   :align: center
   :width: 550px

   Network Model Classes

By convention, network model variables are named ``nm`` and network model class names begin with ``mp.net_model``.


Building a Network Model
^^^^^^^^^^^^^^^^^^^^^^^^

A network model object is created in two steps. The first is to call the constructor of the desired network model class, without arguments. This initializes the :attr:`element_classes <mp.element_container.element_classes>` property with a list of network model element classes. This list can be modified before the second step, which is to call the :meth:`build() <mp.net_model.build>` method, passing in the data model object.

.. _code_net_model_build:
.. code-block::

   nm = mp.net_model_acp();
   nm.build(dm);

The :meth:`build() <mp.net_model.build>` method proceeds through the following stages sequentially, looping through each element at each stage.

   1. **Create** – Instantiate each element object.
   2. **Count and add** - For each element object, determine the number of online elements from the corresponding data model element and, if nonzero, store it in the object and add the object to the :attr:`elements <mp.element_container.elements>` property of the ``nm``.
   3. **Add nodes** – Allow each element to add network nodes, then add voltage variables for each node.
   4. **Add states** – Allow each element to add non-voltage states, then add non-voltage variables for each state.
   5. **Build parameters** – Construct the formulation-specific model parameters for each element, including mappings of element port to network node and element non-voltage state to system non-voltage variable. Add ports to the container object for each element to track per-element port indexing.


Node Types
^^^^^^^^^^

Most problems require that certain nodes be given special treatment depending on their *type*. For example, in the power flow problem, there is typically a single **reference** node, some **PV** nodes, with the rest being **PQ** nodes.

In the current design, each node-creating network model element class implements a :meth:`node_types() <mp.nm_element.node_types>` method that returns information about the types of the nodes it creates. The container object :meth:`node_types() <mp.nm_element.node_types>` method assembles that information for the full set of network nodes. It can also optionally, assign a new reference node if one does not exist. There are also methods, namely :meth:`set_node_type_ref() <mp.nm_element.set_node_type_ref>`, :meth:`set_node_type_pv() <mp.nm_element.set_node_type_pv>`, :meth:`set_node_type_pq() <mp.nm_element.set_node_type_pq>`, for setting the type of a network node and having the relevant elements update their corresponding data model elements.


.. _sec_nm_element:

Network Model Elements
----------------------

A network model element object encapsulates all of the network model parameters for a particular element type. All network model element classes inherit from :class:`mp.nm_element` and also, like the container, from a formulation-specific subclass of :class:`mp.form`. Each element type typically implements its own subclasses, which are further subclassed per formulation. A given network model element object contains the aggregate network model parameters for *all* online instances of that element type, stored in the set of matrices and vectors that correspond to the formulation, e.g. :math:`\BB`, :math:`\KK` and :math:`\pv` from :eq:`eq_GkPx` for DC and :math:`\YY`, :math:`\LL`, :math:`\MM`, :math:`\NN`, :math:`\iv`, and :math:`\sv` from :eq:`eq_Ilin` and :eq:`eq_Slin` for AC.

So, for example, in a system with 1000 in-service transmission lines, the :math:`\YY` parameter in the corresponding AC network model element object would be a 2000 :math:`\times` 2000 matrix for an aggregate 2000-port element, representing the 1000 two-port transmission lines.

By convention, network model element variables are named ``nme`` and network model element class names begin with ``mp.nme``. :numref:`fig_net_model_element_classes` shows the inheritance relationships between a few example network model element classes. Here the :class:`mp.nme_bus_acp` and :class:`mp.nme_gen_acp` classes are used for all problems with an AC polar formulation, while the AC cartesian and DC formulations use their own respective subclasses.

.. _fig_net_model_element_classes:
.. figure:: figures/net-model-element-classes.*
   :alt: Network Model Element Classes
   :align: center

   Network Model Element Classes


Example Elements
^^^^^^^^^^^^^^^^

Here are brief descriptions of the network models for a few simple element types. There are other elements, and the point is that new elements are relatively simple to implement, simply by specifying the nodes, ports and states they add, and the parameters that define the relationships between the states and the port injections.

Bus
'''

A **bus** element inherits from :class:`mp.nme_bus` and defines a single node per in-service bus, with no ports or non-voltage states. So it has no model parameters.


Generator
'''''''''

A **gen** element is a 1-port element that inherits from :class:`mp.nme_gen` and defines a single non-voltage state per in-service generator to represent the power injection. It connects to the node corresponding to a particular bus. The only non-zero parameters are :math:`\KK` (DC) or :math:`\NN` (AC), which are negative identity matrices, since the power injections (into the element) are the negative of the generated power.


Branch
''''''

A **branch** element is a 2-port element that inherits from :class:`mp.nme_branch` with no nodes or non-voltage states. It connects to nodes corresponding to two  particular buses. The only non-zero parameters are :math:`\BB` and :math:`\pv` (DC), or :math:`\YY` (AC).


Load
''''

A **load** element is a 1-port element that inherits from :class:`mp.nme_load` with no ports or states. It connects to the node corresponding to a particular bus. For a simple constant power load, the only non-zero parameters are :math:`\pv` (DC) or :math:`\sv` (AC), equal to the power consumed by the load.


Building Element Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Typically, a network model element builds parameters only for its in-service elements, stacking the corresponding parameters into vectors and matrices, with one row per element of that type. For the DC formulation, these are the three parameters :math:`\BB`, :math:`\KK` and :math:`\pv` from :numref:`sec_nm_formulations_dc`. For the AC formulations they are the six parameters, :math:`\YY`, :math:`\LL`, :math:`\MM`, :math:`\NN`, :math:`\iv`, and :math:`\sv` from :numref:`sec_nm_formulations_ac`.

Take, for example, an AC model with two-port transmission lines modeled by a simple series admittance, where the two ports are labeled with :math:`f` and :math:`t`. For line :math:`i` with series admittance :math:`\cscal{y}^i_s`, we have

.. math::
   :label: eq_single_line_y

   \left[\begin{array}{c}
   \cscal{i}^i_f \\
   \cscal{i}^i_t
   \end{array}\right]
   = \left[\begin{array}{cc}
   \cscal{y}^i_s & -\cscal{y}^i_s \\
   -\cscal{y}^i_s & \cscal{y}^i_s
   \end{array}\right]
   \left[\begin{array}{c}
   \cscal{v}^i_f \\
   \cscal{v}^i_t
   \end{array}\right].

The individual admittance parameters for the :math:`n_k` individual lines are then stacked as follows,

.. math::
   :label: eq_stack_y

   \cmat{Y}_s = \left[\begin{array}{cccc}
        \cscal{y}^1_s & & & \\
        & \cscal{y}^2_s & &\\
         &  & \ddots &  \\
        & & & \cscal{y}^{n_k}_s \\
        \end{array}\right],

to form the admittance matrix parameter :math:`\YY` that we see in :eq:`eq_Ilin` for the corresponding element object.

.. math::
   :label: eq_all_lines_y2

   \YY = \left[\begin{array}{cc}
   \cmat{Y}_s & -\cmat{Y}_s \\
   -\cmat{Y}_s & \cmat{Y}_s
   \end{array}\right]

Stacking the individual port current and voltage variables,

.. math::
   :label: eq_stackiv

   \cvec{i}_f = \left[\begin{array}{c}
        \cscal{i}^1_f \\
        \cscal{i}^2_f \\
        \vdots \\
        \cscal{i}^{n_k}_f \\
        \end{array}\right], \;
   \cvec{i}_t = \left[\begin{array}{c}
        \cscal{i}^1_t \\
        \cscal{i}^2_t \\
        \vdots \\
        \cscal{i}^{n_k}_t \\
        \end{array}\right], \;
   \cvec{v}_f = \left[\begin{array}{c}
        \cscal{v}^1_f \\
        \cscal{v}^2_f \\
        \vdots \\
        \cscal{v}^{n_k}_f \\
        \end{array}\right], \;
   \cvec{v}_t = \left[\begin{array}{c}
        \cscal{v}^1_t \\
        \cscal{v}^2_t \\
        \vdots \\
        \cscal{v}^{n_k}_t \\
        \end{array}\right],

results in the port injection currents from :eq:`eq_GkI` for this aggregate element taking the form

.. math::
   :label: eq_all_lines_y

   \GI(\X) &= \Ilin(\X) = \left[\begin{array}{c}
   \cvec{i}_f \\
   \cvec{i}_t
   \end{array}\right]
   = \YY
   \left[\begin{array}{c}
   \cvec{v}_f \\
   \cvec{v}_t
   \end{array}\right]
   = \YY \V.

When building its parameters, each network model element object also defines an element-node incidence matrix :math:`C` for each of its ports and an element-variable incidence matrix :math:`D` for each non-voltage states. For example, a transmission line element would define two :math:`C` matrices, one mapping branches to their corresponding *from* bus and the other to their corresponding *to* bus.


Aggregation
^^^^^^^^^^^

Since the model parameters are consistent across all network model elements for a given formulation, and the connectivity of the elements is captured in the :math:`C` and :math:`D` incidence matrices for each element type, the network model object can assemble the parameters from all elements into a single aggregate network model characterized by parameters of the same form. This aggregate model can then be used to compute port or node injections from the aggregate system state, as well as any needed derivatives of these injection functions.

For more details on how the aggregation is done, see [TN5]_.

