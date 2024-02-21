.. _howto_gizmo_element:

A Test Gizmo Element Type
=========================

This guide provides another example of the implementation of a new element type. It is not as complete as the one in the :ref:`howto_element` guide, but does illustrate a wider range of network modeling capabilities.

|MATPOWER| includes implementations of numerous standard element types such as buses, generators, loads, transmission lines, transformers, etc. For this discussion, we use a fictitious element type, a **gizmo**, that has no physical equivalent, but illustrates a wide range of modeling possibilities. This is element is used in some of the |MATPOWER| tests to verify custom modeling facilities.

3-Port Gizmo
------------

Let's say that a **gizmo** is a 3 port model with 2 non-voltage states. That is, each gizmo connects to 3 different buses and has 2 additional complex state variables, one essentially a current and the other a power flow. And it includes all six types of components that make up the standard AC network model for any element type, as described in the :ref:`sec_nm_formulations_ac` section of the |MATPOWER-Dev-Manual|.

.. _fig_gizmo_model:
.. figure:: figures/gizmo-model.*
   :alt: AC Model for 3 Port Gizmo
   :align: center
   :width: 400px

   AC Model for 3 Port Gizmo

:numref:`fig_gizmo_model` shows a diagram of the 3-port gizmo, with the various gizmo parameters and the various components summarized in :numref:`tab_gizmo_components`, with the flows through them as functions of the voltages and non-state variables.

.. _tab_gizmo_components:
.. list-table:: Gizmo Components/Parameters
   :widths: 18 82
   :header-rows: 1
   :class: longtable

   * - Parameters
     - Flow Through This Component
   * - :math:`\param{\cscal{i}}_g`
     - complex current is constant |br| :math:`\cscal{i}_1 = \param{\cscal{i}}_g`
   * - :math:`\param{\cscal{s}}_g`
     - complex power is constant |br| :math:`\cscal{s}_2 = \param{\cscal{s}}_g`
   * - :math:`\param{\cscal{y}}_1, \param{\cscal{y}}_2`
     - complex current is proportional to voltage difference (constant impedance) |br| :math:`\cscal{i}_3 = \param{\cscal{y}}_1 (\cscal{v}_1 - \cscal{v}_3)` |br| :math:`\cscal{i}_4 = \param{\cscal{y}}_2 \cscal{v}_2`
   * - :math:`\param{\cscal{m}}_1, \param{\cscal{m}}_2`
     - complex power is proportional to voltage difference |br| :math:`\cscal{s}_5 = \param{\cscal{m}}_1 (\cscal{v}_1 - \cscal{v}_2)` |br| :math:`\cscal{s}_6 = \param{\cscal{m}}_2 \cscal{v}_3`
   * - :math:`\param{\cscal{l}}_g`
     - complex current is proportional to non-voltage state |br| :math:`\cscal{i}_7 = \param{\cscal{l}}_g \cscal{z}_2`
   * - :math:`\param{\cscal{n}}_g`
     - complex power is proportional to non-voltage state |br| :math:`\cscal{s}_8 = \param{\cscal{n}}_g \cscal{z}_1`

Creating a new element type involves defining and implmementing the corresponding classes for each layer (data, network, and mathematical) of modeling, as well as a corresponding data converter element for each data format. Each of these classes has a :meth:`name` method that returns the same value, linking them to one another and together defining the gizmo element type. That is, in each type of element class, you will see the following method defined.

::

        function name = name(obj)
            name = 'gizmo';
        end


Data Model Element
------------------

We begin with the data model element for a **gizmo**. The data for each gizmo consists of the three buses to which its ports are connected, values of the two non-voltage states, and the 8 parameters shown in :numref:`fig_gizmo_model` and :numref:`tab_gizmo_components` above. For the complex values, we specify the real and imaginary parts as separate parameters.

This data is stored in the main data table in the :attr:`tab <mp.dm_element.tab>` property of the data model element object of type :class:`mp.dme_gizmo`. This is a MATLAB |table|_ object, or an :class:`mp_table` object if running in Octave. The names of the columns in :attr:`tab <mp.dm_element.tab>` are shown in :numref:`tab_gizmo_data_model` below. Each row in :attr:`tab <mp.dm_element.tab>` corresponds to an individual gizmo, which means there is a single instance of a gizmo data model element object to hold the data for all gizmos in the system.

.. _tab_gizmo_data_model:
.. list-table:: Gizmo Data Model
   :widths: 18 82
   :header-rows: 1
   :class: longtable

   * - Column Names
     - Description
   * - ``bus_1``, ``bus_2``, ``bus_3``
     - bus numbers for the port 1, 2, and 3 connections, respectively
   * - ``Y1r``, ``Y1i``, ``Y2r``, ``Y2i``
     - real and imaginary parts of parameters :math:`\param{\cscal{y}}_1` and :math:`\param{\cscal{y}}_2`, respectively
   * - ``M1r``, ``M1i``, ``M2r``, ``M2i``
     - real and imaginary parts of parameters :math:`\param{\cscal{m}}_1` and :math:`\param{\cscal{m}}_2`, respectively
   * - ``Lr``, ``Li``
     - real and imaginary parts of parameter :math:`\param{\cscal{l}}_g`
   * - ``Ir``, ``Ii``
     - real and imaginary parts of parameter :math:`\param{\cscal{i}}_g`
   * - ``Nr``, ``Ni``
     - real and imaginary parts of parameter :math:`\param{\cscal{n}}_g`
   * - ``Sr``, ``Si``
     - real and imaginary parts of parameter :math:`\param{\cscal{s}}_g`
   * - ``Zr1``, ``Zi1``, ``Zr2``, ``Zi2``
     - real and imaginary parts of non-voltage state variables :math:`\cscal{z}_1` and :math:`\cscal{z}_2`, respectively


:numref:`code_dme_gizmo` shows the source code for :class:`mp.dme_gizmo`. The first thing to notice is that, as with all data model element classes, it inherits from :class:`mp.dm_element`. Please see the :class:`mp.dm_element` reference documentation for an overview of the functionality provided and for more details on the methods overridden by :class:`mp.dme_gizmo`.

.. _code_dme_gizmo:
.. code-block::
   :linenos:
   :caption: :class:`mp.dme_gizmo`

   classdef dme_gizmo < mp.dm_element
       properties
           bus1        %% bus index vector for port 1
           bus2        %% bus index vector for port 2
           bus3        %% bus index vector for port 3
       end     %% properties

       methods
           function name = name(obj)
               name = 'gizmo';
           end

           function label = label(obj)
               label = 'Test Gizmo';
           end

           function label = labels(obj)
               label = 'Test Gizmos';
           end

           function name = cxn_type(obj)
               name = 'bus';
           end

           function name = cxn_idx_prop(obj)
               name = {'bus1', 'bus2', 'bus3'};
           end

           function names = main_table_var_names(obj)
               names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                   {'bus_1', 'bus_2', 'bus_3', 'Y1r', 'Y1i', 'Y2r', 'Y2i', ...
                   'Lr', 'Li', 'Ir', 'Ii', 'M1r', 'M1i', 'M2r', 'M2i', ...
                   'Nr', 'Ni', 'Sr', 'Si', 'Zr1', 'Zi1', 'Zr2', 'Zi2'});
           end

           function obj = initialize(obj, dm)
               initialize@mp.dm_element(obj, dm);  %% call parent

               %% get bus mapping info
               b2i = dm.elements.bus.ID2i;         %% bus num to idx mapping

               %% set bus index vectors for port connectivity
               obj.bus1 = b2i(obj.tab.bus_1);
               obj.bus2 = b2i(obj.tab.bus_2);
               obj.bus3 = b2i(obj.tab.bus_3);
           end

           function obj = update_status(obj, dm)
               %% get bus status info
               bs = dm.elements.bus.tab.status;        %% bus status

               %% update status of gizmoes connected to isolated/offline buses
               obj.tab.status = obj.tab.status & bs(obj.bus1) & ...
                                                 bs(obj.bus2) & ...
                                                 bs(obj.bus3);

               %% call parent to fill in on/off
               update_status@mp.dm_element(obj, dm);
           end
       end     %% methods
   end         %% classdef


For element types that connect to one or more buses, it is typical to define a property for each port in the data model element class. In our case, there are three properties, :attr:`bus1`, :attr:`bus2`, and :attr:`bus3`, which will hold bus index vectors for ports 1, 2 and 3, respectively. That is ``dme.bus2(k)`` will refer to the index of the bus connected to port 2 of the gizmo defined in row *k* of the data table.

The :meth:`name() <mp.dm_element.name>` method returns ``'gizmo'``, the name used internally for this element type. The :meth:`label() <mp.dm_element.label>` and :meth:`labels() <mp.dm_element.labels>` methods provide strings to use for singular and plural user visible labels to use when displaying gizmo elements.

The :meth:`cxn_type() <mp.dm_element.cxn_type>` and :meth:`cxn_idx_prop() <mp.dm_element.cxn_idx_prop>` methods specify that ``'gizmo'`` objects connect to ``'bus'`` objects and the corresponding bus indices for ports 1, 2, and 3, can be found in properties  :attr:`bus1`, :attr:`bus2`, and :attr:`bus3`, respectively.

The names of the columns in gizmo's main data table are defined by the return value of :meth:`main_table_var_names() <mp.dm_element.main_table_var_names>`. Note that it is important to call the parent method to include the column names common to all data model elements (i.e. ``'uid'``, ``'name'``, ``'status'``, ``'source_uid'``).

The :meth:`initialize() <mp.dm_element.initialize>` method takes advantage of the bus ID to bus index mapping available from the ``'bus'`` data model element object to populate the  :attr:`bus1`, :attr:`bus2`, and :attr:`bus3` properties from the corresponding columns in the main data table.

Finally, :meth:`update_status() <mp.dm_element.update_status>` updates the default online/offline status, which has already been initialized from the ``status`` column of the main data table, to remove from service any gizmo that is connected to an offline bus.

Note that both :meth:`initialize() <mp.dm_element.initialize>` and :meth:`update_status() <mp.dm_element.update_status>` rely on the fact that the corresponding methods have already been called for ``'bus'`` objects before ``'gizmo'`` objects. The order corresponds to their order in :attr:`dm.element_classes <mp.element_container.element_classes>` which is determined by the default defined by the data model class and any |MATPOWER| extensions or options used to modify that default.

The :class:`mp.dme_gizmo` class is also where you would override any of the pretty-printing methods to implement gizmo sections in your pretty-printed output. Until such methods are added to this example, you can look at the data model element classes for other element types for examples (e.g. :class:`mp.dme_bus`, :class:`mp.dme_branch`, :class:`mp.dme_gen`, :class:`mp.dme_load`, etc.)

See |dme_gizmo_m|_ for the complete :class:`mp.dme_gizmo` source.


Data Model Converter Element
----------------------------

*(not yet documented)*


Network Model Element
---------------------

Next we define the **gizmo** network model. The focus will be on the AC model with the assumption that both polar and cartesian voltage formulations should be implemented. Because network models are formulation-specific, we will define a class hierarchy for the network model element.

All Formulations
^^^^^^^^^^^^^^^^

All gizmo network model elements will inherit from :class:`mp.nme_gizmo`, shown in :numref:`code_nme_gizmo`, which in turn inherits from :class:`mp.nm_element`. Please see the :class:`mp.nm_element` reference documentation for an overview of the functionality provided and for more details on the methods overridden by :class:`mp.nme_gizmo` and its subclasses.

.. _code_nme_gizmo:
.. code-block::
   :linenos:
   :caption: :class:`mp.nme_gizmo`

   classdef (Abstract) nme_gizmo < mp.nm_element
       methods
           function name = name(obj)
               name = 'gizmo';
           end

           function np = np(obj)
               np = 3;     %% this is a 3 port element
           end

           function nz = nz(obj)
               nz = 2;     %% 2 (possibly complex) non-voltage states per element
           end
       end     %% methods
   end         %% classdef

Once again, :meth:`name() <mp.nm_element.name>` returns the name used internally for this element type, while the :meth:`np() <mp.nm_element.np>` and :meth:`nz() <mp.nm_element.nz>` methods return the number of ports and non-voltage states, respectively. These are shared by all formulations.


AC Formulations
^^^^^^^^^^^^^^^

Anything specific to all AC formulations is included in the abstract class :class:`mp.nme_gizmo_ac`, shown in :numref:`code_nme_gizmo_ac`, which is a subclass of :class:`mp.nme_gizmo`. Any concrete network model element class that inherits from :class:`mp.nme_gizmo_ac` is also expected to be a subclass of a formulation class that inherits from :class:`mp.form_ac`.

.. _code_nme_gizmo_ac:
.. code-block::
   :linenos:
   :caption: :class:`mp.nme_gizmo_ac`

   classdef (Abstract) nme_gizmo_ac < mp.nme_gizmo% & mp.form_ac
       methods
           function obj = add_zvars(obj, nm, dm, idx)
               tab = obj.data_model_element(dm).tab;
               nk = obj.nk;
               switch idx{:}
                   case 1
                       Zmax = ones(nk, 1);
                       Zr   = tab.Zr1;
                       Zi   = tab.Zi1;
                   case 2
                       Zmax = 2 * ones(nk, 1);
                       Zr   = tab.Zr2;
                       Zi   = tab.Zi2;
               end
               vname_r = sprintf('Zr%d_gizmo', idx{:});
               vname_i = sprintf('Zi%d_gizmo', idx{:});
               nm.add_var('zr', vname_r, nk, Zr, -Zmax, Zmax);
               nm.add_var('zi', vname_i, nk, Zi, -Zmax, Zmax);
           end

           function obj = build_params(obj, nm, dm)
               build_params@mp.nme_gizmo(obj, nm, dm);    %% call parent
               tab = obj.data_model_element(dm).tab;
               nk = obj.nk;

               %% collect parameters from data table
               y1 = tab.Y1r + 1j * tab.Y1i;
               y2 = tab.Y2r + 1j * tab.Y2i;
               ll = tab.Lr + 1j * tab.Li;
               ii = tab.Ir + 1j * tab.Ii;
               m1 = tab.M1r + 1j * tab.M1i;
               m2 = tab.M2r + 1j * tab.M2i;
               nn = tab.Nr + 1j * tab.Ni;
               ss = tab.Sr + 1j * tab.Si;
               zz = zeros(nk, 1);

               %% construct model parameters
               j1 = (1:nk);
               j2 = nk+j1;
               j3 = nk+j2;
               obj.Y = sparse( ...
                   [j1 j1 j1 j2 j2 j2 j3 j3 j3]', ...
                   [j1 j2 j3 j1 j2 j3 j1 j2 j3]', ...
                   [y1; zz; -y1; zz; y2; zz; -y1; zz; y1], 3*nk, 3*nk );
               obj.L = sparse( ...
                   [j1 j1 j2 j2 j3 j3 ]', ...
                   [j1 j2 j1 j2 j1 j2 ]', ...
                   [zz; ll; zz; -ll; zz; zz], 3*nk, 2*nk );
               obj.i = [-ii; ii; zz];
               obj.M = sparse( ...
                   [j1 j1 j1 j2 j2 j2 j3 j3 j3]', ...
                   [j1 j2 j3 j1 j2 j3 j1 j2 j3]', ...
                   [m1; -m1; zz; -m1; m1; zz; zz; zz; m2], 3*nk, 3*nk );
               obj.N = sparse( ...
                   [j1 j1 j2 j2 j3 j3 ]', ...
                   [j1 j2 j1 j2 j1 j2 ]', ...
                   [zz; zz; nn; zz; -nn; zz], 3*nk, 2*nk );
               obj.s = [zz; -ss; ss];
           end
       end     %% methods
   end         %% classdef

The first method defined by :class:`mp.nme_gizmo_ac`, namely :meth:`add_zvars() <mp.nm_element.add_zvars>`, adds variables for the real and imaginary parts of the non-voltage state variables, :math:`\cvec{z}_1` and :math:`\cvec{z}_2`, to the network model, constructing the initial values from the appropriate columns in the data table, and including predefined bounds. We arbitrarily define all gizmos such that their :math:`\cscal{z}` variables, :math:`\cscal{z}_1` and :math:`\cscal{z}_2`, obey :math:`-k \le \Re\{\cscal{z}_k\} \le k` and :math:`-k \le \Im\{\cscal{z}_k\} \le k`. Note that the variable named ``Zr1_gizmo`` is vector containing the real part of :math:`\cscal{z}_1` for all gizmos in the network. Because the voltage variable representation is different for cartesian and polar formulations, the implementation of :meth:`add_vvars() <mp.nm_element.add_vvars>` is deferred to the formulation-specific subclasses below.

The second method, :meth:`build_params() <mp.nm_element.build_params>`, first calls its parent to build the incidence matrices :attr:`C <mp.nm_element.C>` and :attr:`D <mp.nm_element.D>`, then constructs the standard AC model parameters from the data model. The AC model and its parameters are described in :ref:`sec_nm_formulations_ac` in the |MATPOWER-Dev-Manual|.

Recall that, if we omit the arbitrary nonlinear injection components, :math:`\Snln(\X)` or :math:`\Inln(\X)`, the standard AC network model for any element type can be defined in terms of the six parameters in the equations below, namely  :math:`\YY`, :math:`\LL`, :math:`\MM`, :math:`\NN`, :math:`\iv`, and :math:`\sv`. 

.. math::
   :label: eq_Ilin_howto_gizmo

   \Ilin(\X) = \YY \V + \LL \Z + \iv

.. math::
   :label: eq_Slin_howto_gizmo

   \Slin(\X) = \MM \V + \NN \Z + \sv

For a single gizmo, based on :numref:`fig_gizmo_model` and :numref:`tab_gizmo_components`, these parameters would be defined as follows.

.. math::
   :label: eq_i_lin_howto_gizmo

   \YY = \left[\begin{array}{ccc}
       \param{\cscal{y}}_1 & 0 & -\param{\cscal{y}}_1 \\
       0 & \param{\cscal{y}}_2 & 0 \\
       \param{-\cscal{y}}_1 & 0 & \param{\cscal{y}}_1
     \end{array}\right], 
   \LL= \left[\begin{array}{cc}
        0 & \param{\cscal{l}}_g  \\
        0 & -\param{\cscal{l}}_g  \\
        0 & 0
     \end{array}\right],
   \iv = \left[\begin{array}{c}
        -\param{\cscal{i}}_g \\
        \param{\cscal{i}}_g \\
        0
     \end{array}\right]

.. math::
   :label: eq_s_lin_howto_gizmo

   \MM = \left[\begin{array}{ccc}
        \param{\cscal{m}}_1 & -\param{\cscal{m}}_1 & 0 \\
        \param{-\cscal{m}}_1 & \param{\cscal{m}}_1 & 0 \\
        0 & 0 & \param{\cscal{m}}_2
     \end{array}\right],
   \NN = \left[\begin{array}{cc}
        0 & 0 \\
        \param{\cscal{n}}_g & 0 \\
        -\param{\cscal{n}}_g & 0
     \end{array}\right],
   \sv = \left[\begin{array}{c}
        0 \\
        -\param{\cscal{s}}_g \\
        \param{\cscal{s}}_g
     \end{array}\right]

However, :meth:`build_params() <mp.nm_element.build_params>` must build stacked versions of these matrix and vector parameters that include all :math:`n_k` gizmos in the system. For the matrix parameters in :eq:`eq_i_lin_howto_gizmo` and :eq:`eq_s_lin_howto_gizmo`, the stacking is done such that each scalar element is replaced by a corresponding :math:`n_k \times n_k` diagonal matrix. For the vector parameters, each scalar element becomes an :math:`n_k \times 1` vector.


AC Cartesian vs Polar Formulations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once the parameters have been built, all of the differences between the cartesian and polar voltage formulations are handled automatically by inheriting from the appropriate formulation class. For the cartesian voltage formulation, we use :class:`mp.nme_gizmo_acc` which inherits from :class:`mp.nme_gizmo_ac` and :class:`mp.form_acc`.

.. _code_nme_gizmo_acc:
.. code-block::
   :linenos:
   :caption: :class:`mp.nme_gizmo_acc`

   classdef nme_gizmo_acc < mp.nme_gizmo_ac & mp.form_acc
   end

For the polar voltage formulation, we use :class:`mp.nme_gizmo_acp` which inherits from :class:`mp.nme_gizmo_ac` and :class:`mp.form_acp`.

.. _code_nme_gizmo_acp:
.. code-block::
   :linenos:
   :caption: :class:`mp.nme_gizmo_acp`

   classdef nme_gizmo_acp < mp.nme_gizmo_ac & mp.form_acp
   end


Mathematical Model Element
--------------------------

Since the gizmo does not introduce any new costs or gizmo-specific contraints, there is no need for an explicit mathematical model element for gizmos.

This is where you would also put any :meth:`data_model_update() <mp.mm_element.data_model_update>` methods, but the gizmo does not implement any.

.. note::

   The non-voltage state variables are not updated for the power flow, and in the OPF they have the hard-coded limits defined above.


Gizmo Extension
---------------

A |MATPOWER| extension that incorporates this new element can be found in |xt_gizmo_m|_.


.. |dme_gizmo_m| replace:: :file:`lib/t/+mp/dme_gizmo.m`
.. _dme_gizmo_m: https://github.com/MATPOWER/matpower/blob/master/lib/t/%2Bmp/dme_gizmo.m
.. |xt_gizmo_m| replace:: :file:`lib/t/+mp/xt_gizmo.m`
.. _xt_gizmo_m: https://github.com/MATPOWER/matpower/blob/master/lib/t/%2Bmp/xt_gizmo.m
.. |table| replace:: :class:`table`
.. _table: https://www.mathworks.com/help/matlab/ref/table.html
