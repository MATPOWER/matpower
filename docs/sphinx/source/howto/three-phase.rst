.. _howto_three_phase:

How to Run a Three-Phase Power Flow
===================================

This guide describes how to run a power flow on a three-phase, unbalanced model.
|MATPOWER| 8 includes the core functionality required to build out three-phase unbalanced modeling. It also includes a few simple prototype three-phase elements: buses, lines, generators, loads, and transformers and even a bus-link element used to connect three-phase and single-phase parts of a model.

.. note::

    Keep in mind these are **prototype** elements. That means that the data formats used, etc. are not final and there are many additions, enhancements and changes planned for future versions. So don't count on future versions to be backward compatible with these prototypes.

The functionality is accessed through the |MATPOWER| extension :class:`mp.xt_3p` and can be used for power flow (PF), continuation power flow (CPF), and optimal power flow (OPF), by passing the extension as input to the :func:`run_pf`, :func:`run_cpf`, and :func:`run_opf` functions, respectively. 

Some relevant reference documentation can be found in :ref:`ref_xt_3p_classes`, and the test code in |t_run_mp_3p_m|_ includes numerous examples of three-phase cases and how to run them.

Below is an example run.

::

    >> mpopt = mpoption('verbose', 2);
    >> run_pf('t_case3p_a', mpopt, 'mpx', mp.xt_3p)
    
    MATPOWER Version 8.0b1+, 15-Feb-2024
    Power Flow -- AC-polar-power formulation
    
     it    max residual        max ∆x
    ----  --------------  --------------
      0      2.375e+00           -    
      1      4.042e-01       2.287e-01
      2      7.179e-02       7.961e-02
      3      4.570e-03       1.967e-02
      4      2.405e-05       1.400e-03
      5      6.090e-10       6.908e-06
    Newton's method converged in 5 iterations.
    PF successful
    
    PF succeeded in 0.23 seconds (0.22 setup + 0.01 solve)
    ================================================================================
    |     System Summary                                                           |
    ================================================================================
      elements                on     off    total
     --------------------- ------- ------- -------
      3-ph Buses                4       -       4
      3-ph Generators           1       -       1
      3-ph Loads                1       -       1
      3-ph Lines                2       -       2
      3-ph Transformers         1       -       1
    
      Total 3-ph generation               6109.9 kW       4206.5 kVAr
      Total 3-ph load                     5450.0 kW       2442.6 kVAr
      Total 3-ph line loss                 561.5 kW       1173.8 kVAr
      Total 3-ph transformer loss           98.4 kW        590.2 kVAr
    
    ================================================================================
    |     3-ph Bus Data                                                            |
    ================================================================================
      3-ph            Phase A Voltage    Phase B Voltage    Phase C Voltage
     Bus ID   Status   (kV)     (deg)     (kV)     (deg)     (kV)     (deg)
    --------  ------  -------  -------   -------  -------   -------  -------
          1      1     7.1996     0.00    7.1996  -120.00    7.1996   120.00
          2      1     7.1637    -0.14    7.1105  -120.18    7.0821   119.26
          3      1     2.3055    -2.26    2.2547  -123.62    2.2028   114.79
          4      1     2.1750    -4.12    1.9298  -126.80    1.8327   102.85
    
    ================================================================================
    |     3-ph Generator Data                                                      |
    ================================================================================
      3-ph      3-ph             Phase A Power     Phase B Power     Phase C Power
     Gen ID    Bus ID   Status   (kW)    (KVAr)    (kW)    (kVAr)    (kW)    (kVAr)
    --------  --------  ------  -------  ------   -------  ------   -------  ------
          1         1      1    1341.42  970.52   2096.10 1341.41   2672.34 1894.59
    
    ================================================================================
    |     3-ph Load Data                                                           |
    ================================================================================
      3-ph      3-ph             Phase A Power     Phase B Power     Phase C Power
    Load ID    Bus ID   Status   (kW)     (PF)     (kW)     (PF)     (kW)     (PF)
    --------  --------  ------  -------  ------   -------  ------   -------  ------
          1         4      1    1275.00  0.8500   1800.00  0.9000   2375.00  0.9500
    
    ================================================================================
    |     3-ph Line Data                                                           |
    ================================================================================
    -->  Current Injections at "From" Bus
      3-ph    3-ph Bus  3-ph Bus          Phase A Current  Phase B Current  Phase C Current
    Line ID    From ID   To ID    Status   (A)    (deg)     (A)    (deg)     (A)    (deg)
    --------  --------  --------  ------  ------  ------   ------  ------   ------  ------
          1         1         2      1    229.97   -35.9   345.66  -152.6   455.00    84.7
          2         3         4      1    689.64   -35.9  1036.34  -152.6  1364.10    84.7
    
    <--  Current Injections at "To" Bus
      3-ph    3-ph Bus  3-ph Bus          Phase A Current  Phase B Current  Phase C Current
    Line ID    From ID   To ID    Status   (A)    (deg)     (A)    (deg)     (A)    (deg)
    --------  --------  --------  ------  ------  ------   ------  ------   ------  ------
          1         1         2      1    230.06   144.1   345.72    27.4   455.06   -95.3
          2         3         4      1    689.66   144.1  1036.35    27.4  1364.11   -95.3
    
    -->  Power Injections at "From" Bus
      3-ph    3-ph Bus  3-ph Bus          Phase A Power    Phase B Power    Phase C Power
    Line ID    From ID   To ID    Status   (kW)   (kVAr)    (kW)   (kVAr)    (kW)   (kVAr)
    --------  --------  --------  ------  ------  ------   ------  ------   ------  ------
          1         1         2      1    1341.4   970.5   2096.1  1341.4   2672.3  1894.6
          2         3         4      1    1323.5   881.1   2043.4  1133.3   2598.7  1508.6
    
    <--  Power Injections at "To" Bus
      3-ph    3-ph Bus  3-ph Bus          Phase A Power    Phase B Power    Phase C Power
    Line ID    From ID   To ID    Status   (kW)   (kVAr)    (kW)   (kVAr)    (kW)   (kVAr)
    --------  --------  --------  ------  ------  ------   ------  ------   ------  ------
          1         1         2      1   -1337.2  -963.4  -2074.4 -1319.1  -2652.4 -1830.6
          2         3         4      1   -1275.0  -790.2  -1800.0  -871.8  -2375.0  -780.6
    
    ================================================================================
    |     3-ph Transformer Data                                                    |
    ================================================================================
    -->  Current Injections at "From" Bus
      3-ph    3-ph Bus  3-ph Bus          Phase A Current  Phase B Current  Phase C Current
    Xfrm ID    From ID   To ID    Status   (A)    (deg)     (A)    (deg)     (A)    (deg)
    --------  --------  --------  ------  ------  ------   ------  ------   ------  ------
          1         2         3      1    230.06   -35.9   345.72  -152.6   455.06    84.7
    
    <--  Current Injections at "To" Bus
      3-ph    3-ph Bus  3-ph Bus          Phase A Current  Phase B Current  Phase C Current
    Xfrm ID    From ID   To ID    Status   (A)    (deg)     (A)    (deg)     (A)    (deg)
    --------  --------  --------  ------  ------  ------   ------  ------   ------  ------
          1         2         3      1    689.64   144.1  1036.34    27.4  1364.10   -95.3
    
    -->  Power Injections at "From" Bus
      3-ph    3-ph Bus  3-ph Bus          Phase A Power    Phase B Power    Phase C Power
    Xfmr ID    From ID   To ID    Status   (kW)   (kVAr)    (kW)   (kVAr)    (kW)   (kVAr)
    --------  --------  --------  ------  ------  ------   ------  ------   ------  ------
          1         2         3      1    1337.2   963.4   2074.4  1319.1   2652.4  1830.6
    
    <--  Power Injections at "To" Bus
      3-ph    3-ph Bus  3-ph Bus          Phase A Power    Phase B Power    Phase C Power
    Xfmr ID    From ID   To ID    Status   (kW)   (kVAr)    (kW)   (kVAr)    (kW)   (kVAr)
    --------  --------  --------  ------  ------  ------   ------  ------   ------  ------
          1         2         3      1   -1323.5  -881.1  -2043.4 -1133.3  -2598.7 -1508.6


Data Format
-----------

The prototype examples included in |MATPOWER| 8 simply define a few extra fields in the standard |MATPOWER| case struct for the three-phase elements. The example above uses |t_case3p_a|_. This is a case with 4 three-phase buses, connected via 2 lines and 1 transformer, with a three-phase generator representing the substation at bus 1 and a three-phase load at bus 4.

::

    function mpc = t_case3p_a
    %T_CASE3P_A   Four bus, unbalanced 3-phase test case.
    %
    % This data comes from 4Bus-YY-UnB.DSS, a modified version (with unbalanced
    % load) of 4Bus-YY-Bal.DSS [1], the OpenDSS 4 bus IEEE test case with
    % grounded-wye to grounded-wye transformer.
    %
    % [1] https://sourceforge.net/p/electricdss/code/HEAD/tree/trunk/Distrib/IEEETestCases/4Bus-YY-Bal/4Bus-YY-Bal.DSS
    
    %% MATPOWER Case Format : Version 2
    mpc.version = '2';
    
    %%-----  Power Flow Data  -----%%
    %% system MVA base
    mpc.baseMVA = 100;
    
    mpc.bus = [];
    mpc.gen = [];
    mpc.branch = [];
    mpc.gencost = [];
    
    
    %%-----  3 Phase Model Data  -----%%
    %% system data
    mpc.freq = 60;      %% frequency, Hz
    mpc.basekVA = 1000; %% system kVA base
    
    %% bus data
    %	busid	type	basekV	Vm1	Vm2	Vm3	Va1	Va2	Va3
    mpc.bus3p = [
        1	3	12.47	1	1	1	0	-120	120;
        2	1	12.47	1	1	1	0	-120	120;
        3	1	4.16	1	1	1	0	-120	120;
        4	1	4.16	1	1	1	0	-120	120;
    ];
    
    %% branch data
    %	brid	fbus	tbus	status	lcid	len
    mpc.line3p = [
        1	1	2	1	1	2000/5280;
        2	3	4	1	1	2500/5280;
    ];
    
    %% transformer
    %	xfid	fbus	tbus	status	R	X	basekVA	basekV
    mpc.xfmr3p = [
        1	2	3	1	0.01	0.06	6000	12.47;
    ];
    
    %% load
    %	ldid	ldbus	status	Pd1	Pd2	Pd3	ldpf1	ldpf2	ldpf3
    mpc.load3p = [
        1	4	1	1275	1800	2375	0.85	0.9	0.95;
    ];
    
    %% gen
    %	genid	gbus	status	Vg1	Vg2	Vg3	Pg1	Pg2	Pg3	Qg1	Qg2	Qg3
    mpc.gen3p = [
        1	1	1	1	1	1	2000	2000	2000	0	0	0;
    ];
    
    %% line construction
    %	lcid	R11	R21	R31	R22	R32	R33	X11	X21	X31	X22	X32	X33	C11	C21	C31	C22	C32	C33
    mpc.lc = [
        1	0.457541	0.15594 	0.153474	0.466617	0.157996	0.461462	1.078	0.501648	0.384909	1.04813	0.423624	1.06502	15.0671	-4.86241	-1.85323	15.875	-3.09098	14.3254
    ];

The data for the individual elements are found in the fields ``bus3p``, ``line3p``, ``xfrm3p``, ``load3p``, and ``gen3p`` fields of ``mpc``. The ``lc`` field contains the per-mile impedance parameters for different line construction configurations referenced by the individual lines.

Besides the fields for the individual elements, cases with three-phase elements are also expected to include ``freq``, the system frequency in Hertz, and ``base_kVA``, the system per-unit kVA base for three-phase portions of the network.

Details of the data for the individual three-phase elements are summarized in the tables below.


``bus3p``
^^^^^^^^^

===  =========  =====================
Col  Name       Description — *see also* :class:`mp.dme_bus3p`
===  =========  =====================
  1  busid      unique 3-phase bus ID (positive integer)
  2  type       bus type: 1 = PQ, 2 = PV, 3 = reference (voltage ref + slack)
  3  basekV     nominal voltage in kV
  4  Vm1        phase A voltage magnitude in p.u.
  5  Vm2        phase B voltage magnitude in p.u.
  6  Vm3        phase C voltage magnitude in p.u.
  7  Va1        phase A voltage angle in degrees (nominal 0)
  8  Va2        phase B voltage angle in degrees (nominal -120)
  9  Va3        phase C voltage angle in degrees (nominal 120)
===  =========  =====================


``line3p``
^^^^^^^^^^

===  =========  =====================
Col  Name       Description — *see also* :class:`mp.dme_line3p`
===  =========  =====================
  1  brid       unique 3-phase line ID (positive integer)
  2  fbus       "from" bus ID
  3  tbus       "to" bus ID
  4  status     1 = in-service, 0 = out-of-service
  5  lcid       line construction ID (to look up in ``lc`` table)
  6  len        line length in miles
===  =========  =====================


``xfmr3p``
^^^^^^^^^^

===  =========  =====================
Col  Name       Description — *see also* :class:`mp.dme_xfmr3p`
===  =========  =====================
  1  xfid       unique 3-phase transformer ID (positive integer)
  2  fbus       "from" bus ID
  3  tbus       "to" bus ID
  4  status     1 = in-service, 0 = out-of-service
  5  R          resistance (p.u. on transformer base)
  6  X          reactance (p.u. on transformer base)
  7  basekVA    transformer per-unit power base in kVA
  8  basekV     transformer per-unit voltage base in kV
===  =========  =====================


``load3p``
^^^^^^^^^^

===  =========  =====================
Col  Name       Description — *see also* :class:`mp.dme_load3p`
===  =========  =====================
  1  ldid       unique 3-phase load ID (positive integer)
  2  ldbus      bus ID
  3  status     1 = in-service, 0 = out-of-service
  4  Pd1        phase A active power demand in kW
  5  Pd2        phase B active power demand in kW
  6  Pd3        phase C active power demand in kW
  7  lpf1       phase A load power factor
  8  lpf2       phase B load power factor
  9  lpf3       phase C load power factor
===  =========  =====================


``gen3p``
^^^^^^^^^

===  =========  =====================
Col  Name       Description — *see also* :class:`mp.dme_gen3p`
===  =========  =====================
  1  genid      unique 3-phase generator ID (positive integer)
  2  gbus       bus ID
  3  status     1 = in-service, 0 = out-of-service
  4  Vg1        phase A voltage setpoint in p.u.
  5  Vg2        phase B voltage setpoint in p.u.
  6  Vg3        phase C voltage setpoint in p.u.
  7  Pg1        phase A active power injection in kW
  8  Pg2        phase B active power injection in kW
  9  Pg3        phase C active power injection in kW
 10  Qg1        phase A reactive power injection in kVAr
 11  Qg2        phase B reactive power injection in kVAr
 12  Qg3        phase C reactive power injection in kVAr
===  =========  =====================


``lc``
^^^^^^

===  =========  =====================
Col  Name       Description — *see also* :class:`mp.dme_line3p`
===  =========  =====================
  1  lcid       unique 3-phase line construction record ID (positive integer)
  2  R11        resistance (p.u. per mile), element (1,1) of 3x3 matrix
  3  R21        resistance (p.u. per mile), elements (2,1) & (1,2) of 3x3 matrix
  4  R31        resistance (p.u. per mile), elements (3,1) & (1,3) of 3x3 matrix
  5  R22        resistance (p.u. per mile), element (2,2) of 3x3 matrix
  6  R32        resistance (p.u. per mile), elements (3,2) & (2,3) of 3x3 matrix
  7  R33        resistance (p.u. per mile), element (3,3) of 3x3 matrix
  8  X11        reactance (p.u. per mile), element (1,1) of 3x3 matrix)
  9  X21        reactance (p.u. per mile), elements (2,1) & (1,2) of 3x3 matrix
 10  X31        reactance (p.u. per mile), elements (3,1) & (1,3) of 3x3 matrix
 11  X22        reactance (p.u. per mile), element (2,2) of 3x3 matrix
 12  X32        reactance (p.u. per mile), elements (3,2) & (2,3) of 3x3 matrix
 13  X33        reactance (p.u. per mile), element (3,3) of 3x3 matrix
 14  C11        capacitance (nF per mile), element (1,1) of 3x3 matrix
 15  C21        capacitance (nF per mile), elements (2,1) & (1,2) of 3x3 matrix
 16  C31        capacitance (nF per mile), elements (3,1) & (1,3) of 3x3 matrix
 17  C22        capacitance (nF per mile), element (2,2) of 3x3 matrix
 18  C32        capacitance (nF per mile), elements (3,2) & (2,3) of 3x3 matrix
 19  C33        capacitance (nF per mile), element (3,3) of 3x3 matrix
===  =========  =====================


``buslink``
^^^^^^^^^^^

===  =========  =====================
Col  Name       Description — *see also* :class:`mp.dme_buslink`
===  =========  =====================
  1  linid      unique 3-phase buslink ID (positive integer)
  2  busid      bus ID of single-phase bus
  3  bus3pid    bus ID of three-phase bus
  4  status     1 = in-service, 0 = out-of-service
===  =========  =====================


Example Cases
-------------

The following cases are included in |lib_t|_.

============= ================== ======= ======= ================
Name          Description        3-phase 1-phase bus links
============= ================== ======= ======= ================
|t_case3p_a|_ 4 bus 3-phase case 4       --      --
|t_case3p_b|_ 6 bus hybrid case  4       2       1 PQ
|t_case3p_c|_ 6 bus hybrid case  4       2       1 PV (1-phase side)
|t_case3p_d|_ 6 bus hybrid case  4       2       1 PV (3-phase side)
|t_case3p_e|_ 5 bus hybrid case  4       1       1 REF (1-phase side)
|t_case3p_f|_ 21 bus hybrid case 12      9       3 PQ
|t_case3p_g|_ 21 bus hybrid case 12      9       3 REF-PQ, PV-PQ, PQ-PQ
|t_case3p_h|_ 21 bus hybrid case 12      9       3 REF-PQ, PQ-PV, PQ-PQ
============= ================== ======= ======= ================

The data for the 4 bus, three-phase system in |t_case3p_a|_ comes from ``4Bus-YY-UnB.DSS``, a modified version (with unbalanced load) of |4Bus_YY_Bal_DSS|_, the OpenDSS 4-bus IEEE test case with grounded-wye to grounded-wye transformer. [#]_

The five and six bus cases connect this 4-bus three-phase case to 1 or 2 single-phase buses. And the 21-bus cases are based on connecting three copies of the 4 bus three-phase case to the single-phase |case9|_.

.. |t_run_mp_3p_m| replace:: :file:`t_run_mp_3p.m`
.. |lib_t| replace:: :file:`lib/t`
.. |case9| replace:: :file:`case9.m`
.. |4Bus_YY_Bal_DSS| replace:: :file:`4Bus-YY-Bal.DSS`
.. |t_case3p_a| replace:: :file:`t_case3p_a.m`
.. |t_case3p_b| replace:: :file:`t_case3p_b.m`
.. |t_case3p_c| replace:: :file:`t_case3p_c.m`
.. |t_case3p_d| replace:: :file:`t_case3p_d.m`
.. |t_case3p_e| replace:: :file:`t_case3p_e.m`
.. |t_case3p_f| replace:: :file:`t_case3p_f.m`
.. |t_case3p_g| replace:: :file:`t_case3p_g.m`
.. |t_case3p_h| replace:: :file:`t_case3p_h.m`

.. [#] https://sourceforge.net/p/electricdss/code/HEAD/tree/trunk/Distrib/IEEETestCases/4Bus-YY-Bal/4Bus-YY-Bal.DSS
.. _t_run_mp_3p_m: https://github.com/MATPOWER/matpower/blob/master/lib/t/t_run_mp_3p.m
.. _lib_t: https://github.com/MATPOWER/matpower/tree/master/lib/t
.. _case9: https://github.com/MATPOWER/matpower/blob/master/data/case9.m
.. _4Bus_YY_Bal_DSS: https://sourceforge.net/p/electricdss/code/HEAD/tree/trunk/Distrib/IEEETestCases/4Bus-YY-Bal/4Bus-YY-Bal.DSS
.. _t_case3p_a: https://github.com/MATPOWER/matpower/blob/master/lib/t/t_case3p_a.m
.. _t_case3p_b: https://github.com/MATPOWER/matpower/blob/master/lib/t/t_case3p_b.m
.. _t_case3p_c: https://github.com/MATPOWER/matpower/blob/master/lib/t/t_case3p_c.m
.. _t_case3p_d: https://github.com/MATPOWER/matpower/blob/master/lib/t/t_case3p_d.m
.. _t_case3p_e: https://github.com/MATPOWER/matpower/blob/master/lib/t/t_case3p_e.m
.. _t_case3p_f: https://github.com/MATPOWER/matpower/blob/master/lib/t/t_case3p_f.m
.. _t_case3p_g: https://github.com/MATPOWER/matpower/blob/master/lib/t/t_case3p_g.m
.. _t_case3p_h: https://github.com/MATPOWER/matpower/blob/master/lib/t/t_case3p_h.m
