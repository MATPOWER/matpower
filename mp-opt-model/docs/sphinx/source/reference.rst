.. _sec_mpom_reference:

Reference
=========


|MPOM| Classes
--------------

Use :class:`mp.opt_model` to build and solve your mathematical programming and optimization models.

.. toctree::

   classes/mp/opt_model
   classes/mp/set_manager
   classes/mp/set_manager_opt_model
   classes/mp/sm_variable
   classes/mp/sm_lin_constraint
   classes/mp/sm_nln_constraint
   classes/mp/sm_nln_cost
   classes/mp/sm_quad_cost
   classes/mp/sm_quad_constraint

Deprecated:

.. toctree::

   classes/opt_model
   classes/mp_idx_manager
   classes/mp/sm_quad_cost_legacy


|MPOM| Functions
----------------

.. _mpom_master_functions:

Master Functions
^^^^^^^^^^^^^^^^

The |MPOM| master functions provide unified interfaces to multiple solvers for linear program (LP), quadratic program (QP), quadratically-constrained quadratic program (QCQP), nonlinear program (NLP), mixed-integer program (MILP/MIQP) optimization, as well as nonlinear equation (NLEQ) solving and parameterized nonlinear equation (PNE) solution tracing.

.. toctree::

   functions/miqps_master
   functions/nleqs_master
   functions/nlps_master
   functions/pnes_master
   functions/qcqps_master
   functions/qps_master


Utility Functions
^^^^^^^^^^^^^^^^^

Use these functions to convert constraints or to copy data from one struct to another.

.. toctree::

   functions/convert_constraint_multipliers
   functions/convert_lin_constraint_multipliers
   functions/convert_lin_constraint
   functions/convert_quad_constraint
   functions/nested_struct_copy
   functions/mp/struct2object


Options Handling Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^

Use these functions to set up input options for individual solvers.

.. toctree::

   functions/artelys_knitro_options
   functions/cplex_options
   functions/glpk_options
   functions/gurobi_options
   functions/highs_options
   functions/ipopt_options
   functions/mosek_options
   functions/mosek_symbcon
   functions/osqp_options
   functions/mpopt2nleqopt
   functions/mpopt2nlpopt
   functions/mpopt2pneopt
   functions/mpopt2qcqpopt
   functions/mpopt2qpopt


Version Information Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use these functions to check the version of the installed |MPOM|, Gurobi, HiGHS, and OSQP.

.. toctree::

   functions/gurobiver
   functions/highsver
   functions/knitrover
   functions/mpomver
   functions/osqpver


Solver Interfaces & Implementations
-----------------------------------

These functions provide the implementations and interfaces for the solvers accessible through the :ref:`mpom_master_functions` above.

MIQP Solver
^^^^^^^^^^^

.. toctree::

   functions/miqps_cplex
   functions/miqps_glpk
   functions/miqps_gurobi
   functions/miqps_highs
   functions/miqps_mosek
   functions/miqps_ot


NLEQ Solver
^^^^^^^^^^^

.. toctree::

   functions/nleqs_core
   functions/nleqs_fd_newton
   functions/nleqs_fsolve
   functions/nleqs_gauss_seidel
   functions/nleqs_newton


NLP Solver
^^^^^^^^^^

.. toctree::

   functions/nlp_consfcn
   functions/nlp_costfcn
   functions/nlp_hessfcn
   functions/nlps_fmincon
   functions/nlps_ipopt
   functions/nlps_knitro


PNE Solver
^^^^^^^^^^

.. toctree::

   functions/pne_callback_default
   functions/pne_callback_nose
   functions/pne_callback_target_lam
   functions/pne_detect_events
   functions/pne_detected_event
   functions/pne_event_nose
   functions/pne_event_target_lam
   functions/pne_pfcn_arc_len
   functions/pne_pfcn_natural
   functions/pne_pfcn_pseudo_arc_len
   functions/pne_register_callbacks
   functions/pne_register_events


LP/QP Solver
^^^^^^^^^^^^

.. toctree::

   functions/qps_bpmpd
   functions/qps_clp
   functions/qps_cplex
   functions/qps_glpk
   functions/qps_gurobi
   functions/qps_highs
   functions/qps_ipopt
   functions/qps_knitro
   functions/qps_mosek
   functions/qps_osqp
   functions/qps_ot


QCQP Solver
^^^^^^^^^^^

.. toctree::

   functions/qcqp_nlp_consfcn
   functions/qcqp_nlp_costfcn
   functions/qcqp_nlp_hessfcn
   functions/qcqps_gurobi
   functions/qcqps_knitro
   functions/qcqps_nlps


|MPOM| Examples
---------------

These are examples of using |MPOM| to build and solve various types of problems.

.. toctree::

   functions/lp_ex1
   functions/milp_ex1
   functions/nleqs_master_ex1
   functions/nleqs_master_ex2
   functions/nlps_master_ex1
   functions/nlps_master_ex2
   functions/pne_ex1
   functions/qcqp_ex1
   functions/qp_ex1

The following Live Scripts, also found in the ``lib/t`` directory, help to
illustrate MP-Opt-Model usage.

- `milp_example1.mlx`_ - *Good introduction to build and solve a model with MP-Opt-Model.*
- `qcqp_example1.mlx`_ - *Example of quadratic constraints and QCQP.*

.. _milp_example1.mlx: https://github.com/MATPOWER/mp-opt-model/blob/master/lib/t/milp_example1.mlx
.. _qcqp_example1.mlx: https://github.com/MATPOWER/mp-opt-model/blob/master/lib/t/qcqp_example1.mlx


|MPOM| Tests
------------

These functions test that |MPOM| is installed and functioning as expected.

.. toctree::

   functions/test_mp_opt_model
   functions/t_nested_struct_copy
   functions/t_nleqs_master
   functions/t_pnes_master
   functions/t_qps_master
   functions/t_qcqps_master
   functions/t_miqps_master
   functions/t_nlps_master
   functions/t_mp_opt_model
   functions/t_mm_solve_leqs
   functions/t_mm_solve_nleqs
   functions/t_mm_solve_pne
   functions/t_mm_solve_qps
   functions/t_mm_solve_qcqps
   functions/t_mm_solve_miqps
   functions/t_mm_solve_nlps

Legacy Tests
^^^^^^^^^^^^

Tests of deprecated functionality.

.. toctree::

   functions/t_have_fcn
   functions/t_opt_model
   functions/t_om_solve_leqs
   functions/t_om_solve_nleqs
   functions/t_om_solve_pne
   functions/t_om_solve_qps
   functions/t_om_solve_qcqps
   functions/t_om_solve_miqps
   functions/t_om_solve_nlps


Private Functions
-----------------

The following are private functions that implement detection of specific
optional functionality. They are not intended to be called directly, but
rather are used to extend the capabilities of :func:`have_feature`.

.. toctree::

   functions/have_feature_bpmpd
   functions/have_feature_catchme
   functions/have_feature_clp
   functions/have_feature_cplex
   functions/have_feature_evalc
   functions/have_feature_fmincon_ipm
   functions/have_feature_fmincon
   functions/have_feature_fsolve
   functions/have_feature_glpk
   functions/have_feature_gurobi
   functions/have_feature_highs
   functions/have_feature_intlinprog
   functions/have_feature_ipopt_auxdata
   functions/have_feature_ipopt
   functions/have_feature_isequaln
   functions/have_feature_knitro
   functions/have_feature_knitromatlab
   functions/have_feature_linprog_ds
   functions/have_feature_linprog
   functions/have_feature_mosek
   functions/have_feature_opti_clp
   functions/have_feature_optim
   functions/have_feature_optimoptions
   functions/have_feature_osqp
   functions/have_feature_quadprog_ls
   functions/have_feature_quadprog
   functions/have_feature_sdpt3
   functions/have_feature_sedumi
   functions/have_feature_yalmip
