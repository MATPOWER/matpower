Classes
=======

Task Classes
------------

Core Task Classes
+++++++++++++++++

.. toctree::
   :name: sec_task_classes

   mp/task
   mp/task_pf
   mp/task_cpf
   mp/task_opf

Legacy Task Classes
+++++++++++++++++++

Used by MP-Core when called by the *legacy* |/MATPOWER/| *framework*.

.. toctree::
   :name: sec_legacy_task_classes

   mp/task_pf_legacy
   mp/task_cpf_legacy
   mp/task_opf_legacy
   mp/task_shared_legacy


Data Model Classes
------------------

Containers
++++++++++

.. toctree::
   :name: sec_data_model_container_classes

   mp/data_model
   mp/data_model_cpf
   mp/data_model_opf

Elements
++++++++

.. toctree::
   :name: sec_data_model_element_classes

   mp/dm_element
   mp/dme_branch
   mp/dme_branch_opf
   mp/dme_bus
   mp/dme_bus_opf
   mp/dme_gen
   mp/dme_gen_opf
   mp/dme_load
   mp/dme_load_cpf
   mp/dme_load_opf
   mp/dme_shunt_cpf
   mp/dme_shunt
   mp/dme_shunt_opf


Element Mixins
++++++++++++++

.. toctree::
   :name: sec_data_model_element_mixin_classes

   mp/dme_shared_opf


Data Model Converter Classes
----------------------------

Containers
++++++++++

.. toctree::
   :name: sec_data_model_converter_container_classes

   mp/dm_converter
   mp/dm_converter_mpc2
   mp/dm_converter_mpc2_legacy

Elements
++++++++

.. toctree::
   :name: sec_data_model_converter_element_classes

   mp/dmc_element
   mp/dmce_branch_mpc2
   mp/dmce_bus_mpc2
   mp/dmce_gen_mpc2
   mp/dmce_load_mpc2
   mp/dmce_shunt_mpc2


Network Model Classes
---------------------

Containers
++++++++++

.. toctree::
   :name: sec_net_model_container_classes

   mp/form
   mp/form_ac
   mp/form_acc
   mp/form_acp
   mp/form_dc
   mp/net_model
   mp/net_model_ac
   mp/net_model_acc
   mp/net_model_acp
   mp/net_model_dc

Elements
++++++++

.. toctree::
   :name: sec_net_model_element_classes

   mp/nm_element
   mp/nme_branch
   mp/nme_branch_ac
   mp/nme_branch_acc
   mp/nme_branch_acp
   mp/nme_branch_dc
   mp/nme_bus
   mp/nme_bus_acc
   mp/nme_bus_acp
   mp/nme_bus_dc
   mp/nme_gen
   mp/nme_gen_ac
   mp/nme_gen_acc
   mp/nme_gen_acp
   mp/nme_gen_dc
   mp/nme_load
   mp/nme_load_ac
   mp/nme_load_acc
   mp/nme_load_acp
   mp/nme_load_dc
   mp/nme_shunt
   mp/nme_shunt_ac
   mp/nme_shunt_acc
   mp/nme_shunt_acp
   mp/nme_shunt_dc


Mathematical Model Classes
--------------------------

Containers
++++++++++

.. toctree::
   :name: sec_math_model_container_classes

   mp/math_model
   mp/math_model_pf
   mp/math_model_pf_ac
   mp/math_model_pf_acci
   mp/math_model_pf_accs
   mp/math_model_pf_acpi
   mp/math_model_pf_acps
   mp/math_model_pf_dc
   mp/math_model_cpf_acc
   mp/math_model_cpf_acci
   mp/math_model_cpf_accs
   mp/math_model_cpf_acp
   mp/math_model_cpf_acpi
   mp/math_model_cpf_acps
   mp/math_model_opf
   mp/math_model_opf_ac
   mp/math_model_opf_acc
   mp/math_model_opf_acci
   mp/math_model_opf_acci_legacy
   mp/math_model_opf_accs
   mp/math_model_opf_accs_legacy
   mp/math_model_opf_acp
   mp/math_model_opf_acpi
   mp/math_model_opf_acpi_legacy
   mp/math_model_opf_acps
   mp/math_model_opf_acps_legacy
   mp/math_model_opf_dc
   mp/math_model_opf_dc_legacy

Container Mixins
++++++++++++++++

.. toctree::
   :name: sec_math_model_mixin_classes

   mp/mm_shared_pfcpf
   mp/mm_shared_pfcpf_ac
   mp/mm_shared_pfcpf_ac_i
   mp/mm_shared_pfcpf_acc
   mp/mm_shared_pfcpf_acci
   mp/mm_shared_pfcpf_accs
   mp/mm_shared_pfcpf_acp
   mp/mm_shared_pfcpf_acpi
   mp/mm_shared_pfcpf_acps
   mp/mm_shared_pfcpf_dc
   mp/mm_shared_opf_legacy

Elements
++++++++

.. toctree::
   :name: sec_math_model_element_classes

   mp/mm_element
   mp/mme_branch
   mp/mme_branch_pf_ac
   mp/mme_branch_pf_dc
   mp/mme_branch_opf
   mp/mme_branch_opf_ac
   mp/mme_branch_opf_acc
   mp/mme_branch_opf_acp
   mp/mme_branch_opf_dc
   mp/mme_bus
   mp/mme_bus_pf_ac
   mp/mme_bus_pf_dc
   mp/mme_bus_opf_ac
   mp/mme_bus_opf_acc
   mp/mme_bus_opf_acp
   mp/mme_bus_opf_dc
   mp/mme_gen
   mp/mme_gen_pf_ac
   mp/mme_gen_pf_dc
   mp/mme_gen_opf
   mp/mme_gen_opf_ac
   mp/mme_gen_opf_dc
   mp/mme_load
   mp/mme_load_pf_ac
   mp/mme_load_pf_dc
   mp/mme_load_cpf
   mp/mme_shunt
   mp/mme_shunt_pf_ac
   mp/mme_shunt_pf_dc
   mp/mme_shunt_cpf


Miscellaneous Classes
---------------------

.. toctree::
   :name: sec_misc_classes

   mp_table
   mp_table_subclass
   mp/cost_table
   mp/cost_table_utils
   mp/element_container
   mp/mapped_array
   mp/NODE_TYPE


|MATPOWER| Extension Classes
----------------------------

Base
++++

.. toctree::
   :name: sec_xt_classes

   mp/extension

.. _ref_xt_reserves_classes:

OPF Fixed Zonal Reserves Extension
++++++++++++++++++++++++++++++++++

.. toctree::
   :name: sec_xt_reserves_class

   mp/xt_reserves

Other classes belonging to :class:`mp.xt_reserves` extension:
   .. toctree::
      :name: sec_xt_reserves_classes

      mp/dmce_reserve_gen_mpc2
      mp/dmce_reserve_zone_mpc2
      mp/dme_reserve_gen
      mp/dme_reserve_zone
      mp/mme_reserve_gen
      mp/mme_reserve_zone

.. _ref_xt_3p_classes:

Three-Phase Prototype Extension
+++++++++++++++++++++++++++++++

.. toctree::
   :name: sec_xt_3p_class

   mp/xt_3p

Data model converter element classes belonging to :class:`mp.xt_3p` extension:
   .. toctree::
      :name: sec_xt_3p_dmce_classes

      mp/dmce_bus3p_mpc2
      mp/dmce_gen3p_mpc2
      mp/dmce_load3p_mpc2
      mp/dmce_line3p_mpc2
      mp/dmce_xfmr3p_mpc2
      mp/dmce_buslink_mpc2

Data model element classes belonging to :class:`mp.xt_3p` extension:
   .. toctree::
      :name: sec_xt_3p_dme_classes

      mp/dme_bus3p
      mp/dme_gen3p
      mp/dme_load3p
      mp/dme_line3p
      mp/dme_xfmr3p
      mp/dme_buslink
      mp/dme_bus3p_opf
      mp/dme_gen3p_opf
      mp/dme_load3p_opf
      mp/dme_line3p_opf
      mp/dme_xfmr3p_opf
      mp/dme_buslink_opf

Network model element classes belonging to :class:`mp.xt_3p` extension:
   .. toctree::
      :name: sec_xt_3p_nme_classes

      mp/nme_bus3p
      mp/nme_bus3p_acc
      mp/nme_bus3p_acp
      mp/nme_gen3p
      mp/nme_gen3p_acc
      mp/nme_gen3p_acp
      mp/nme_load3p
      mp/nme_line3p
      mp/nme_xfmr3p
      mp/nme_buslink
      mp/nme_buslink_acc
      mp/nme_buslink_acp

Mathematical model element classes belonging to :class:`mp.xt_3p` extension:
   .. toctree::
      :name: sec_xt_3p_mme_classes

      mp/mme_bus3p
      mp/mme_gen3p
      mp/mme_line3p
      mp/mme_xfmr3p
      mp/mme_buslink
      mp/mme_buslink_pf_ac
      mp/mme_buslink_pf_acc
      mp/mme_buslink_pf_acp
      mp/mme_bus3p_opf_acc
      mp/mme_bus3p_opf_acp
      mp/mme_gen3p_opf
      mp/mme_line3p_opf
      mp/mme_xfmr3p_opf
      mp/mme_buslink_opf
      mp/mme_buslink_opf_acc
      mp/mme_buslink_opf_acp

.. _ref_xt_legacy_dcline_classes:

Legacy DC Line Extension
++++++++++++++++++++++++

For more details, see :ref:`howto_element`.

.. toctree::
   :name: sec_xt_legacy_dcline_class

   mp/xt_legacy_dcline

Data model converter element class belonging to :class:`mp.xt_legacy_dcline` extension:
   .. toctree::
      :name: sec_xt_legacy_dcline_dmce_classes

      mp/dmce_legacy_dcline_mpc2

Data model element classes belonging to :class:`mp.xt_legacy_dcline` extension:
   .. toctree::
      :name: sec_xt_legacy_dcline_dme_classes

      mp/dme_legacy_dcline
      mp/dme_legacy_dcline_opf

Network model element classes belonging to :class:`mp.xt_legacy_dcline` extension:
   .. toctree::
      :name: sec_xt_legacy_dcline_nme_classes

      mp/nme_legacy_dcline
      mp/nme_legacy_dcline_ac
      mp/nme_legacy_dcline_acc
      mp/nme_legacy_dcline_acp
      mp/nme_legacy_dcline_dc

Mathematical model element classes belonging to :class:`mp.xt_legacy_dcline` extension:
   .. toctree::
      :name: sec_xt_legacy_dcline_mme_classes

      mp/mme_legacy_dcline
      mp/mme_legacy_dcline_pf_ac
      mp/mme_legacy_dcline_pf_dc
      mp/mme_legacy_dcline_opf
      mp/mme_legacy_dcline_opf_ac
      mp/mme_legacy_dcline_opf_dc

.. _ref_xt_oval_cap_curve_classes:

Example User Constraint Extension
+++++++++++++++++++++++++++++++++

For more details, see :ref:`howto_add_constraint`.

.. toctree::
   :name: sec_xt_oval_cap_curve

   mp/xt_oval_cap_curve

Mathematical model element class belonging to :class:`mp.xt_oval_cap_curve` extension:
   .. toctree::
      :name: xt_oval_cap_curve_mme_classes

      mp/mme_gen_opf_ac_oval
