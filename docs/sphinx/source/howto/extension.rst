.. _howto_extension:

How to Create an Extension
==========================

This guide demonstrates how to create a |MATPOWER| extension to package up your new functionality and make it accessible to a |MATPOWER| end user.

In general, new or modified functionality takes the form of additional classes that inherit from or override functionality in existing classes, or classes that are entirely new. The :ref:`Customizing MATPOWER <sec_customizing>` section in the |MATPOWER-Dev-Manual| describes how |MATPOWER| selects the classes to use for a given run and some mechanisms to make changes to the defaults. In particular, to change the element classes associated with a particular container (e.g. data model, network model, etc.) |MATPOWER| uses *element class modifiers* which are summarized in the :ref:`sec_element_classes` section and the :ref:`tab_element_class_modifiers` table.

In this guide, we will use the example of adding a new network element, namely the model of a DC line used by legacy |MATPOWER| versions. The classes that implement this functionality are described in :ref:`howto_element`. Here we show how to package these into a |MATPOWER| extension.

This extension adds a new ``legacy_dcline`` element for power flow (PF), continuation power flow (CPF), and optimal power flow (OPF) problems. For PF and OPF, it includes both DC and AC formulations, and for the AC formulations both cartesian and polar voltage representations are implemented. :numref:`code_xt_legacy_dcline` shows the code for this extension. As with all |MATPOWER| extensions, it inherits from :class:`mp.extension`. In this example, no changes are required for the task or container classes, so only the methods with names ending in ``_element_classes`` are overridden.

For each model layer, it adds a single element class. For the data model element, that class depends on the task. For the network model element, it depends on the formulation, and for the math model element, it depends on both.


.. _code_xt_legacy_dcline:
.. code-block::
   :linenos:
   :caption: :class:`mp.xt_legacy_dcline`

   classdef xt_legacy_dcline < mp.extension
       methods
           function dmc_elements = dmc_element_classes(obj, dmc_class, fmt, mpopt)
               switch fmt
                   case 'mpc2'
                       dmc_elements = { @mp.dmce_legacy_dcline_mpc2 };
                   otherwise
                       dmc_elements = {};
               end
           end
   
           function dm_elements = dm_element_classes(obj, dm_class, task_tag, mpopt)
               switch task_tag
                   case {'PF', 'CPF'}
                       dm_elements = { @mp.dme_legacy_dcline };
                   case 'OPF'
                       dm_elements = { @mp.dme_legacy_dcline_opf };
                   otherwise
                       dm_elements = {};
               end
           end
   
           function nm_elements = nm_element_classes(obj, nm_class, task_tag, mpopt)
               switch task_tag
                   case {'PF', 'CPF'}
                       v_cartesian = mpopt.pf.v_cartesian;
                   case {'OPF'}
                       v_cartesian = mpopt.opf.v_cartesian;
               end
               switch upper(mpopt.model)
                   case 'AC'
                       if v_cartesian
                           nm_elements = { @mp.nme_legacy_dcline_acc };
                       else
                           nm_elements = { @mp.nme_legacy_dcline_acp };
                       end
                   case 'DC'
                       nm_elements = { @mp.nme_legacy_dcline_dc };
                   otherwise
                       nm_elements = {};
               end
           end
   
           function mm_elements = mm_element_classes(obj, mm_class, task_tag, mpopt)
               switch task_tag
                   case {'PF', 'CPF'}
                       switch upper(mpopt.model)
                           case 'AC'
                               mm_elements = { @mp.mme_legacy_dcline_pf_ac };
                           case 'DC'
                               mm_elements = { @mp.mme_legacy_dcline_pf_dc };
                       end
                   case {'OPF'}
                       switch upper(mpopt.model)
                           case 'AC'
                               mm_elements = { @mp.mme_legacy_dcline_opf_ac };
                           case 'DC'
                               mm_elements = { @mp.mme_legacy_dcline_opf_dc };
                       end
                   otherwise
                       dm_elements = {};
               end
           end
       end     %% methods
   end         %% classdef

See |xt_legacy_dcline_m|_ for the complete, documented :class:`mp.dmce_legacy_dcline_mpc2` source.


Using the Extension
-------------------

To make use of the extension, simply pass it, as an optional argument, immediately following ``'mpx'`` to :func:`run_pf`, :func:`run_cpf`, or :func:`run_opf`, as follows.

::

    >> mpopt = mpoption('verbose', 0);
    >> run_opf('t_case9_dcline', mpopt, 'mpx', mp.xt_legacy_dcline)
    
    OPF succeeded in 0.27 seconds (0.23 setup + 0.04 solve)
    Objective Function Value = 6518.89 $/hr
    ================================================================================
    |     System Summary                                                           |
    ================================================================================
      elements                on     off    total
     --------------------- ------- ------- -------
      Buses                     9       -       9
        Areas                                   1
        Zones                                   1
      Generators                3       -       3
      Loads                     3       -       3
      Branches                  9       -       9
        Lines                   9       -       9
        Transformers            0       -       0
      DC Lines                  3       1       4
    
      Total generation                     319.4 MW          1.5 MVAr
      Total max generation capacity        820.0 MW        900.0 MVAr
      Total min generation capacity        110.0 MW       -900.0 MVAr
      Total load                           315.0 MW        115.0 MVAr
      Total branch losses                    3.00 MW      -124.10 MVAr
      Total DC line losses                   1.40 MW       -10.59 MVAr
    
                                               minimum                        maximum
                                   -----------------------------  -----------------------------
      Bus voltage magnitude              1.066 p.u. @ bus 30            1.100 p.u. @ bus 1
      Bus voltage angle                   -4.51 deg @ bus 5               4.11 deg @ bus 2
      Lambda P (LMP active power)       14.99 $/MWh @ bus 6            24.48 $/MWh @ bus 9
      Lambda Q (LMP reactive power)   -0.62 $/MVArh @ bus 30          0.43 $/MVArh @ bus 7
    
    ================================================================================
    |     Bus Data                                                                 |
    ================================================================================
                          Voltage            Lambda (LMP)
     Bus ID   Status  Mag(pu)  Ang(deg)  P($/MWh)  Q($/MVAr-hr)
    --------  ------  -------  --------  --------  ------------
          1      1     1.100     0.000    15.954         0.000
          2      1     1.100     4.107    24.035         0.000
         30      1     1.066     2.277    15.047        -0.623
          4      1     1.094    -2.470    15.967         0.298
          5      1     1.078    -4.508    15.952        -0.000
          6      1     1.085    -0.277    14.992        -0.600
          7      1     1.081    -2.160    24.476         0.427
          8      1     1.097     0.205    24.043         0.112
          9      1     1.067    -4.470    24.476        -0.079
    
    ================================================================================
    |     Load Data                                                                |
    ================================================================================
                                 Power Consumption
    Load ID    Bus ID   Status   P (MW)   Q (MVAr)
    --------  --------  ------  --------  --------
          1         5      1       90.0      30.0
          2         7      1      100.0      35.0
          3         9      1      125.0      50.0
                                --------  --------
                       Total:     315.0     115.0
    
    ================================================================================
    |     Branch Data                                                              |
    ================================================================================
     Branch     From       To             From Bus Injection   To Bus Injection
       ID      Bus ID    Bus ID   Status   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)
    --------  --------  --------  ------  --------  --------  --------  --------
          1         1         4      1      90.00     14.19    -90.00    -10.24
          2         4         5      1      47.54      1.31    -47.20    -18.11
          3         5         6      1     -48.85    -12.48     49.68    -25.81
          4        30         6      1      88.02    -32.65    -88.02     37.19
          5         6         7      1      38.35    -11.39    -38.20    -11.87
          6         7         8      1     -69.64    -23.13     70.01      8.58
          7         8         2      1    -131.39     -0.99    131.39      9.95
          8         8         9      1      61.37     -7.59    -60.34    -23.03
          9         9         4      1     -51.06    -36.97     51.36     18.92
    
    ================================================================================
    |     DC Line Data                                                             |
    ================================================================================
     DC Line    From       To              Power Flow (MW)      Loss    Reactive Inj (MVAr)
       ID      Bus ID    Bus ID   Status    From       To       (MW)      From       To
    --------  --------  --------  ------  --------  --------  --------  --------  --------
          1        30         4      1      10.00      8.90      1.10    -10.00     10.00
          2         7         9      1       7.84      7.84      0.00     -0.00     -0.00
          3         5         8      0       0.00      0.00      0.00      0.00      0.00
          4         5         9      1       6.06      5.75      0.30     -0.59    -10.00
    
    ================================================================================
    |     Bus Constraints                                                          |
    ================================================================================
                         Voltage Magnitude Limits
     Bus ID     mu LB      LB       vm       UB       mu UB
    --------  ---------  -------  -------  -------   --------
          1       -       0.900    1.100    1.100    566.231
          2       -       0.900    1.100    1.100    197.010
    
    ================================================================================
    |     Generator Constraints                                                    |
    ================================================================================
                                      Active Power Limits
     Gen ID    Bus ID     mu LB      LB       pg       UB       mu UB
    --------  --------  ---------  -------  -------  -------   --------
          1         1      9.046    90.00    90.00   250.00       -   
          2        30       -       10.00    98.02   270.00      0.047
    
                                     Reactive Power Limits
     Gen ID    Bus ID     mu LB      LB       qg       UB       mu UB
    --------  --------  ---------  -------  -------  -------   --------
          2        30      0.623  -300.00   -22.65   300.00       -   
    
    ================================================================================
    |     Branch Constraints            (S in MVA)                                 |
    ================================================================================
     Branch     From        "From" End       Limit       "To" End         To
       ID      Bus ID    mu_sm_fr   sm_fr    sm_ub    sm_to    mu_sm_to  Bus ID
    --------  --------  ---------  -------  -------  -------  ---------  --------
          5         6      2.762    40.00    40.00    40.00      7.375         7
    
    ================================================================================
    |     DC Line Constraints                                                      |
    ================================================================================
     DC Line    From       To                   Active Power Flow (MW)
       ID      Bus ID    Bus ID     mu LB       LB      p_fr      UB      mu UB
    --------  --------  --------  ---------  -------  -------  -------  ---------
          1        30         4       -        1.00    10.00    10.00      0.760


See Also
--------

- :ref:`howto_element`
- :ref:`sec_customizing` in the |MATPOWER-Dev-Manual|
- :ref:`sec_extensions` in the |MATPOWER-Dev-Manual|
- :ref:`ref_xt_legacy_dcline_classes` (:class:`mp.xt_legacy_dcline`)
- :ref:`ref_xt_reserves_classes` (:class:`mp.xt_reserves`)
- :ref:`ref_xt_3p_classes` (:class:`mp.xt_3p`)


.. |xt_legacy_dcline_m| replace:: :file:`lib/t/+mp/xt_legacy_dcline.m`
.. _xt_legacy_dcline_m: https://github.com/MATPOWER/matpower/blob/master/lib/t/%2Bmp/xt_legacy_dcline.m
