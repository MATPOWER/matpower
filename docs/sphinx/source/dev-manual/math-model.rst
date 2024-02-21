.. _sec_math_model:

Mathematical Model Object
=========================

The mathematical model, or math model, formulates and defines the mathematical problem to be solved. That is, it determines the variables, constraints, and objective that define the problem. This takes on different forms depending on the task and the formulation.

Power Flow
    The *power flow* problem involves solving a system of nonlinear equations for the vector :math:`\x`.

    .. math::
       :label: eq_pf_ac_form

       \rvec{f}(\x) = \rvec{0}

    For the DC version, the function :math:`\rvec{f}(\x)` is linear, so the problem takes the more specific form,

    .. math::
       :label: eq_pf_dc_form

       \param{\rmat{A}} \x - \param{\rvec{b}} = \rvec{0}.

Continuation Power Flow
    The *continuation power flow* problem involves tracing the solution curve for a parameterized system of equations, as the parameter :math:`\lambda` is varied.

    .. math::
       :label: eq_cpf_form

       \rvec{f}(\x, \lambda) = \rvec{0},

Optimal Power Flow
    The *optimal power flow* problem, on the other hand, is a constrained  optimization problem of the form,

    .. math::
       :label: eq_opf_ac_form

       \min_\x \rvec{f}(\x) \\
       \textrm{such that} \qquad \rvec{g}(\x) &= \rvec{0} \\
       \rvec{h}(\x) &\le \rvec{0} \\
       \param{\x}_{\min} \le \x &\le \param{\x}_{\max}.

    This reduces to a simple quadratic program (QP) for the DC OPF case,

    .. math::
       :label: eq_opf_dc_form

       \min_\x \trans{\x} \param{\rmat{Q}} \x &+ \trans{\param{\rvec{c}}} \x + \param{k} \\
       \textrm{such that} \qquad \param{\rvec{l}} &\le \param{\rmat{A}} \x \le \param{\rvec{u}} \\
       \param{\x}_{\min} & \le \x \le \param{\x}_{\max}.


Mathematical Models
-------------------

A math model object is a container for math model element objects and it is also an |MPOM>| object. All math model classes inherit from :class:`mp.math_model` and therefore also from :class:`mp.element_container`, :class:`opt_model`, and  :class:`mp_idx_manager`. Concrete math model classes are task and formulation specific as illustrated in :numref:`fig_math_model_classes`, and sometimes inherit from abstract mix-in classes that are shared across tasks or formulations. These shared classes are described further in :numref:`sec_math_model_shared`.

.. _fig_math_model_classes:
.. figure:: figures/math-model-classes.*
   :alt: Math Model Classes
   :align: center

   Math Model Classes

By convention, math model variables are named ``mm`` and math model class names begin with ``mp.math_model``.


Building a Mathematical Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A math model object is created in two steps. The first is to call the constructor of the desired math model class, without arguments. This initializes the :attr:`element_classes <mp.element_container.element_classes>` property with a list of math model element classes. This list can be modified before the second step, which is to call the :meth:`build() <mp.math_model.build>` method, passing in the network and data model objects and a |MATPOWER| options struct.

.. _code_math_model_build:
.. code-block::

   mm = mp.math_model_opf_acps();
   mm.build(nm, dm, mpopt);

The :meth:`build() <mp.math_model.build>` method proceeds through the following stages sequentially, looping through each element for the last 3 stages.

   1. **Create** – Instantiate each element object.
   2. **Count and add** - For each element object, determine the number of online elements from the corresponding data model element and, if nonzero, add the object to the :attr:`elements <mp.element_container.elements>` property of the ``mm``.
   3. **Add auxiliary data** – Add auxiliary data, e.g. network node types, for use by the model.
   4. **Add variables** – Add variables and allow each element to add their own variables to the model.
   5. **Add constraints** – Add constraints and allow each element to add their own constraints to the model.
   6. **Add costs** – Add costs and allow each element to add their own costs to the model.

The adding of variables, constraints and costs to the model is done by the math model and model model element objects using the interfaces provided by |MPOM>|.


Solving a Math Model
^^^^^^^^^^^^^^^^^^^^

Once the math model build is complete and it contains the full set of variables, constraints and costs for the model, the solver options are initialized by calling the :meth:`solve_opts() <mp.math_model.solve_opts>` method and then passed to the :meth:`solve` method.

.. _code_math_model_solve:
.. code-block::

   opt = mm.solve_opts(nm, dm, mpopt);
   mm.solve(opt);

The :meth:`solve` method, also inherited from |MPOM>|, invokes the appropriate solver based on the characteristics of the model and the options provided.


Updating Network and Data Models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The solved math model can then be used to update the solved state of the network and data models by calling the :meth:`network_model_x_soln() <mp.math_model.network_model_x_soln>` and :meth:`data_model_update() <mp.math_model.data_model_update>` methods, respectively.

.. _code_math_model_updates:
.. code-block::

   nm = mm.network_model_x_soln(nm);
   dm = mm.data_model_update(nm, dm, mpopt);

The math model's :meth:`data_model_update() <mp.math_model.data_model_update>` method cycles through the math model element objects, calling the :meth:`data_model_update() <mp.mm_element.data_model_update>` for each element.


.. _sec_mm_element:

Mathematical Model Elements
---------------------------

A math model element object typically does not contain any data, but only the methods that are used to build the math model and update the corresponding data model element once the math model has been solved.

All math model element classes inherit from :class:`mp.mm_element`. Each element type typically implements its own subclasses, which are further subclassed where necessary per task and formulation, as with the container class.

By convention, math model element variables are named ``mme`` and math model element class names begin with ``mp.mme``. :numref:`fig_math_model_element_classes` shows the inheritance relationships between a few example math model element classes. Here the :class:`mp.nme_bus_pf_acp` and :class:`mp.nme_bus_opf_acp` classes are used for PF and OPF problems, respectively, with an AC polar formulation. AC cartesian and DC formulations use their own respective task-specific subclasses. And each element type, has a similar set of task and formulation-specific subclasses, such as those for :class:`mp.mme_gen`.

.. _fig_math_model_element_classes:
.. figure:: figures/math-model-element-classes.*
   :alt: Math Model Element Classes
   :align: center
   :width: 650px

   Math Model Element Classes


Adding Variables, Constraints, and Costs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Both the ``mm`` container object and the ``mme`` element objects can add their own variables, costs and constraints to the model.

For a standard optimal power flow, for example, the optimization variables are added by the container object, since they are determined directly from state variables of the *(container)* network model object. Similarly, the nodal power or current balance constraints are added by the container since they are built directly from the port injection functions of the aggregate network model.

However, generator cost functions and any variables and constraints associated with piecewise linear generator costs are added by the appropriate subclass of :class:`mp.mme_gen`, since they relate only to generator model parameters. Similarly, branch flow and branch angle difference constraints are added by the appropriate subclass of :class:`mp.mme_branch`, since they are specific to branches and are completely independent of other element types.


Updating Data Model Elements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The data in the data model is stored primarily in its individual element objects, so it makes sense that the individual math model element objects would be responsible for extracting the math model solution data relevant to a given element and updating the corresponding data model element. This updating is performed by the :meth:`data_model_update() <mp.mm_element.data_model_update>` method.

The updating of each data model element is done in two steps. First :meth:`data_model_update() <mp.mm_element.data_model_update>` calls :meth:`data_model_update_off() <mp.mm_element.data_model_update_off>` to handle any offline units (e.g. to zero out any solution values), then :meth:`data_model_update_on() <mp.mm_element.data_model_update_on>` to handle the online units.

For example, updating the branch power flows and shadow prices on the flow and angle difference limits in the branch data model element is done by :meth:`data_model_update_on() <mp.mm_element.data_model_update_on>` in the appropriate subclass of :class:`mp.mme_branch`.


.. _sec_math_model_shared:


Shared Classes
--------------

In some cases, there is code shared between math model classes across differnt tasks, e.g. PF and CPF. In order to avoid code duplication, another hierarchy of abstract mix-in classes is used to implement methods for this shared functionality. By convention, the names of these classes begin with ``mp.mm_shared_``.

For example, a method to evaluate the node balance equations and corresponding Jacobian are used by both the PF and CPF. Putting this method in a shared class, allows its functionality to be inherited by concrete math model classes for both PF and CPF.
