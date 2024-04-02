Legacy
======

This section contains reference documentation for the **legacy** |*MATPOWER*| **framework** (see :ref:`sec_two_frameworks` in the |MATPOWER-Dev-Manual|) and the rest of the legacy codebase inherited from |MATPOWER| 7 and earlier.

Legacy Class
------------

.. toctree::

   classes/opf_model


Legacy Functions
----------------

Top-Level Simulation Functions
++++++++++++++++++++++++++++++

.. toctree::

   functions/runpf
   functions/runcpf
   functions/runopf
   functions/runuopf
   functions/rundcpf
   functions/rundcopf
   functions/runduopf
   functions/runopf_w_res


Input/Output Functions
++++++++++++++++++++++

.. toctree::

   functions/caseformat
   functions/cdf2mpc
   functions/loadcase
   functions/mpoption
   functions/printpf
   functions/psse2mpc
   functions/save2psse
   functions/savecase
   functions/savechgtab


Data Conversion Functions
+++++++++++++++++++++++++

.. toctree::

   functions/ext2int
   functions/e2i_data
   functions/e2i_field
   functions/int2ext
   functions/i2e_data
   functions/i2e_field
   functions/get_reorder
   functions/set_reorder


Power Flow Functions
++++++++++++++++++++

.. toctree::

   functions/calc_v_i_sum
   functions/calc_v_pq_sum
   functions/calc_v_y_sum
   functions/dcpf
   functions/fdpf
   functions/gausspf
   functions/make_vcorr
   functions/make_zpv
   functions/newtonpf
   functions/newtonpf_I_cart
   functions/newtonpf_I_hybrid
   functions/newtonpf_I_polar
   functions/newtonpf_S_cart
   functions/newtonpf_S_hybrid
   functions/order_radial
   functions/pfsoln
   functions/radial_pf
   functions/zgausspf


Continuation Power Flow Functions
+++++++++++++++++++++++++++++++++

.. toctree::

   functions/cpf_corrector
   functions/cpf_current_mpc
   functions/cpf_default_callback
   functions/cpf_detect_events
   functions/cpf_flim_event
   functions/cpf_flim_event_cb
   functions/cpf_nose_event
   functions/cpf_nose_event_cb
   functions/cpf_p
   functions/cpf_p_jac
   functions/cpf_plim_event
   functions/cpf_plim_event_cb
   functions/cpf_predictor
   functions/cpf_qlim_event
   functions/cpf_qlim_event_cb
   functions/cpf_register_callback
   functions/cpf_register_event
   functions/cpf_tangent
   functions/cpf_target_lam_event
   functions/cpf_target_lam_event_cb
   functions/cpf_vlim_event
   functions/cpf_vlim_event_cb


OPF and Wrapper Functions
+++++++++++++++++++++++++

.. toctree::

   functions/opf
   functions/dcopf
   functions/fmincopf
   functions/uopf


Other OPF Functions
+++++++++++++++++++

.. toctree::

   functions/dcopf_solver
   functions/nlpopf_solver
   functions/makeAang
   functions/makeApq
   functions/makeAvl
   functions/makeAy
   functions/margcost
   functions/opf_args
   functions/opf_setup
   functions/opf_execute
   functions/opf_branch_ang_fcn
   functions/opf_branch_ang_hess
   functions/opf_branch_flow_fcn
   functions/opf_branch_flow_hess
   functions/opf_current_balance_fcn
   functions/opf_current_balance_hess
   functions/opf_gen_cost_fcn
   functions/opf_legacy_user_cost_fcn
   functions/opf_power_balance_fcn
   functions/opf_power_balance_hess
   functions/opf_veq_fcn
   functions/opf_veq_hess
   functions/opf_vlim_fcn
   functions/opf_vlim_hess
   functions/opf_vref_fcn
   functions/opf_vref_hess
   functions/totcost
   functions/update_mupq


OPF User Callback Functions
+++++++++++++++++++++++++++

.. toctree::

   functions/add_userfcn
   functions/remove_userfcn
   functions/run_userfcn
   functions/toggle_dcline
   functions/toggle_iflims
   functions/toggle_reserves
   functions/toggle_softlims


Power Flow Derivative Functions
+++++++++++++++++++++++++++++++

.. toctree::

   functions/dIbr_dV
   functions/dSbr_dV
   functions/dAbr_dV
   functions/dImis_dV
   functions/dSbus_dV
   functions/d2Ibr_dV2
   functions/d2Sbr_dV2
   functions/d2Abr_dV2
   functions/d2Imis_dV2
   functions/d2Imis_dVdSg
   functions/d2Sbus_dV2


LP, QP, MILP & MIQP Solver Functions
++++++++++++++++++++++++++++++++++++

.. toctree::

   functions/miqps_matpower
   functions/qps_matpower


Matrix Building Functions
+++++++++++++++++++++++++

.. toctree::

   functions/makeB
   functions/makeBdc
   functions/makeJac
   functions/makeLODF
   functions/makePTDF
   functions/makeSbus
   functions/makeSdzip
   functions/makeYbus


Utility Functions
+++++++++++++++++

.. toctree::

   functions/apply_changes
   functions/bustypes
   functions/calc_branch_angle
   functions/case_info
   functions/compare_case
   functions/define_constants
   functions/extract_islands
   functions/feval_w_path
   functions/find_bridges
   functions/find_islands
   functions/genfuels
   functions/gentypes
   functions/get_losses
   functions/hasPQcap
   functions/idx_brch
   functions/idx_bus
   functions/idx_cost
   functions/idx_ct
   functions/idx_dcline
   functions/idx_gen
   functions/isload
   functions/load2disp
   functions/loadshed
   functions/modcost
   functions/mpver
   functions/poly2pwl
   functions/polycost
   functions/pqcost
   functions/scale_load
   functions/total_load


Private Feature Detection Functions
+++++++++++++++++++++++++++++++++++

.. toctree::

   functions/have_feature_e4st
   functions/have_feature_minopf
   functions/have_feature_most
   functions/have_feature_mp_core
   functions/have_feature_pdipmopf
   functions/have_feature_regexp_split
   functions/have_feature_scpdipmopf
   functions/have_feature_sdp_pf
   functions/have_feature_smartmarket
   functions/have_feature_syngrid
   functions/have_feature_table
   functions/have_feature_tralmopf


Other Functions
+++++++++++++++

.. toctree::

   functions/connected_components
   functions/mpoption_info_clp
   functions/mpoption_info_cplex
   functions/mpoption_info_fmincon
   functions/mpoption_info_glpk
   functions/mpoption_info_gurobi
   functions/mpoption_info_intlinprog
   functions/mpoption_info_ipopt
   functions/mpoption_info_knitro
   functions/mpoption_info_linprog
   functions/mpoption_info_mips
   functions/mpoption_info_mosek
   functions/mpoption_info_osqp
   functions/mpoption_info_quadprog
   functions/mpoption_old
   functions/psse_convert
   functions/psse_convert_hvdc
   functions/psse_convert_xfmr
   functions/psse_parse
   functions/psse_parse_line
   functions/psse_parse_section
   functions/psse_read


Legacy Tests
------------

Legacy |MATPOWER| Tests
+++++++++++++++++++++++

.. toctree::

   functions/t_apply_changes
   functions/t_auction_minopf
   functions/t_auction_mips
   functions/t_auction_tspopf_pdipm
   functions/t_chgtab
   functions/t_cpf
   functions/t_dcline
   functions/t_ext2int2ext
   functions/t_feval_w_path
   functions/t_get_losses
   functions/t_hasPQcap
   functions/t_hessian
   functions/t_islands
   functions/t_jacobian
   functions/t_load2disp
   functions/t_loadcase
   functions/t_makeLODF
   functions/t_makePTDF
   functions/t_margcost
   functions/t_miqps_matpower
   functions/t_modcost
   functions/t_mpoption
   functions/t_mpoption_ov
   functions/t_off2case
   functions/t_opf_dc_bpmpd
   functions/t_opf_dc_clp
   functions/t_opf_dc_cplex
   functions/t_opf_dc_default
   functions/t_opf_dc_glpk
   functions/t_opf_dc_gurobi
   functions/t_opf_dc_ipopt
   functions/t_opf_dc_mips
   functions/t_opf_dc_mips_sc
   functions/t_opf_dc_mosek
   functions/t_opf_dc_osqp
   functions/t_opf_dc_ot
   functions/t_opf_default
   functions/t_opf_fmincon
   functions/t_opf_ipopt
   functions/t_opf_knitro
   functions/t_opf_minopf
   functions/t_opf_mips
   functions/t_opf_model
   functions/t_opf_softlims
   functions/t_opf_tspopf_pdipm
   functions/t_opf_tspopf_scpdipm
   functions/t_opf_tspopf_tralm
   functions/t_opf_userfcns
   functions/t_pf_ac
   functions/t_pf_dc
   functions/t_pf_radial
   functions/t_printpf
   functions/t_psse
   functions/t_qps_matpower
   functions/t_runmarket
   functions/t_runopf_w_res
   functions/t_scale_load
   functions/t_total_load
   functions/t_totcost
   functions/t_vdep_load

Legacy |MATPOWER| Test Data
+++++++++++++++++++++++++++

.. toctree::

   functions/opf_nle_fcn1
   functions/opf_nle_hess1
   functions/t_auction_case
   functions/t_case30_userfcns
   functions/t_case9_dcline
   functions/t_case9_opf
   functions/t_case9_opfv2
   functions/t_case9_pf
   functions/t_case9_pfv2
   functions/t_case9_save2psse
   functions/t_case_ext
   functions/t_case_int
   functions/t_cpf_cb1
   functions/t_cpf_cb2

