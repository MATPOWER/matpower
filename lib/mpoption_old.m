function [options, names] = mpoption_old(varargin)
%MPOPTION_OLD  Used to set and retrieve old-style MATPOWER options vector.
%
%   OPT = MPOPTION_OLD
%       returns the default options vector
%
%   OPT = MPOPTION_OLD(NAME1, VALUE1, NAME2, VALUE2, ...)
%       returns the default options vector with new values for up to 7
%       options, NAME# is the name of an option, and VALUE# is the new
%       value.
%
%   OPT = MPOPTION_OLD(OPT, NAME1, VALUE1, NAME2, VALUE2, ...)
%       same as above except it uses the options vector OPT as a base
%       instead of the default options vector.
%
%   Examples:
%       opt = mpoption_old('PF_ALG', 2, 'PF_TOL', 1e-4);
%       opt = mpoption_old(opt, 'OPF_ALG', 565, 'VERBOSE', 2);
%
%   The currently defined options are as follows:
%
%      idx - NAME, default          description [options]
%      ---   -------------          -----------------------------------------
%   power flow options
%       1  - PF_ALG, 1              AC power flow algorithm
%           [   1 - Newton's method                                         ]
%           [   2 - Fast-Decoupled (XB version)                             ]
%           [   3 - Fast-Decoupled (BX version)                             ]
%           [   4 - Gauss-Seidel                                            ]
%       2  - PF_TOL, 1e-8           termination tolerance on per unit
%                                   P & Q mismatch
%       3  - PF_MAX_IT, 10          maximum number of iterations for
%                                   Newton's method
%       4  - PF_MAX_IT_FD, 30       maximum number of iterations for 
%                                   fast decoupled method
%       5  - PF_MAX_IT_GS, 1000     maximum number of iterations for 
%                                   Gauss-Seidel method
%       6  - ENFORCE_Q_LIMS, 0      enforce gen reactive power limits
%                                   at expense of |V|
%           [    0 - do NOT enforce limits                                  ]
%           [    1 - enforce limits, simultaneous bus type conversion       ]
%           [    2 - enforce limits, one-at-a-time bus type conversion      ]
%       10 - PF_DC, 0               DC modeling for power flow & OPF
%           [    0 - use AC formulation & corresponding algorithm options   ]
%           [    1 - use DC formulation, ignore AC algorithm options        ]
%   OPF options
%       11 - OPF_ALG, 0             solver to use for AC OPF
%           [    0 - choose default solver based on availability in the     ]
%           [        following order, 540, 560                              ]
%           [  500 - MINOPF, MINOS-based solver, requires optional          ]
%           [        MEX-based MINOPF package, available from:              ]
%           [        http://www.pserc.cornell.edu/minopf/                   ]
%           [  520 - fmincon, MATLAB Optimization Toolbox >= 2.x            ]
%           [  540 - PDIPM, primal/dual interior point method, requires     ]
%           [        optional MEX-based TSPOPF package, available from:     ]
%           [        http://www.pserc.cornell.edu/tspopf/                   ]
%           [  545 - SC-PDIPM, step-controlled variant of PDIPM, requires   ]
%           [        TSPOPF (see 540)                                       ]
%           [  550 - TRALM, trust region based augmented Langrangian        ]
%           [        method, requires TSPOPF (see 540)                      ]
%           [  560 - MIPS, MATPOWER Interior Point Solver                   ]
%           [        primal/dual interior point method (pure MATLAB)        ]
%           [  565 - MIPS-sc, step-controlled variant of MIPS               ]
%           [        primal/dual interior point method (pure MATLAB)        ]
%           [  580 - IPOPT, requires MEX interface to IPOPT solver          ]
%           [        available from: https://projects.coin-or.org/Ipopt/    ]
%           [  600 - KNITRO, requires MATLAB Optimization Toolbox and       ]
%           [        KNITRO libraries available from: http://www.ziena.com/ ]
%       16 - OPF_VIOLATION, 5e-6    constraint violation tolerance
%       17 - CONSTR_TOL_X, 1e-4     termination tol on x for fmincon/Knitro
%       18 - CONSTR_TOL_F, 1e-4     termination tol on f for fmincon/Knitro
%       19 - CONSTR_MAX_IT, 0       max number of iterations for fmincon
%                                   [       0 => default                    ]
%       24 - OPF_FLOW_LIM, 0        qty to limit for branch flow constraints
%           [   0 - apparent power flow (limit in MVA)                      ]
%           [   1 - active power flow (limit in MW)                         ]
%           [   2 - current magnitude (limit in MVA at 1 p.u. voltage)      ]
%       25 - OPF_IGNORE_ANG_LIM, 0  ignore angle difference limits for branches
%                                   even if specified           [   0 or 1  ]
%       26 - OPF_ALG_DC, 0          solver to use for DC OPF
%           [    0 - choose default solver based on availability in the     ]
%           [        following order: 500, 600, 700, 100, 300, 200          ]
%           [  100 - BPMPD, requires optional MEX-based BPMPD_MEX package   ]
%           [        available from: http://www.pserc.cornell.edu/bpmpd/    ]
%           [  200 - MIPS, MATLAB Interior Point Solver                     ]
%           [        primal/dual interior point method (pure MATLAB)        ]
%           [  250 - MIPS-sc, step-controlled variant of MIPS               ]
%           [  300 - MATLAB Optimization Toolbox, QUADPROG, LINPROG         ]
%           [  400 - IPOPT, requires MEX interface to IPOPT solver          ]
%           [        available from: https://projects.coin-or.org/Ipopt/    ]
%           [  500 - CPLEX, requires MATLAB interface to CPLEX solver       ]
%           [  600 - MOSEK, requires MATLAB interface to MOSEK solver       ]
%           [        available from: http://www.mosek.com/                  ]
%           [  700 - GUROBI, requires Gurobi optimizer (v. 5+)              ]
%           [        available from: http://www.gurobi.com/                 ]
%   output options
%       31 - VERBOSE, 1             amount of progress info printed
%           [   0 - print no progress info                                  ]
%           [   1 - print a little progress info                            ]
%           [   2 - print a lot of progress info                            ]
%           [   3 - print all progress info                                 ]
%       32 - OUT_ALL, -1            controls pretty-printing of results
%           [  -1 - individual flags control what prints                    ]
%           [   0 - do not print anything                                   ]
%           [       (overrides individual flags)                            ]
%           [   1 - print everything                                        ]
%           [       (overrides individual flags)                            ]
%       33 - OUT_SYS_SUM, 1         print system summary        [   0 or 1  ]
%       34 - OUT_AREA_SUM, 0        print area summaries        [   0 or 1  ]
%       35 - OUT_BUS, 1             print bus detail            [   0 or 1  ]
%       36 - OUT_BRANCH, 1          print branch detail         [   0 or 1  ]
%       37 - OUT_GEN, 0             print generator detail      [   0 or 1  ]
%                                   (OUT_BUS also includes gen info)
%       38 - OUT_ALL_LIM, -1        controls what constraint info is printed
%           [  -1 - individual flags control what constraint info prints    ]
%           [   0 - no constraint info (overrides individual flags)         ]
%           [   1 - binding constraint info (overrides individual flags)    ]
%           [   2 - all constraint info (overrides individual flags)        ]
%       39 - OUT_V_LIM, 1           control output of voltage limit info
%           [   0 - do not print                                            ]
%           [   1 - print binding constraints only                          ]
%           [   2 - print all constraints                                   ]
%           [   (same options for OUT_LINE_LIM, OUT_PG_LIM, OUT_QG_LIM)     ]
%       40 - OUT_LINE_LIM, 1        control output of line flow limit info
%       41 - OUT_PG_LIM, 1          control output of gen P limit info
%       42 - OUT_QG_LIM, 1          control output of gen Q limit info
%       44 - OUT_FORCE, 0           print results even if success = 0
%                                                               [   0 or 1  ]
%       52 - RETURN_RAW_DER, 0      return constraint and derivative info
%                                   in results.raw (in fields g, dg, df, d2f)
%   FMINCON options
%       55 - FMC_ALG, 4             algorithm used by fmincon for OPF
%                                   for Optimization Toolbox 4 and later
%            [  1 - active-set                                              ]
%            [  2 - interior-point, w/default 'bfgs' Hessian approx         ]
%            [  3 - interior-point, w/ 'lbfgs' Hessian approx               ]
%            [  4 - interior-point, w/exact user-supplied Hessian           ]
%            [  5 - interior-point, w/Hessian via finite differences        ]
%
%   KNITRO options
%       58 - KNITRO_OPT, 0          a non-zero integer N indicates that all
%                                   KNITRO options should be handled by a
%                                   KNITRO options file named
%                                   'knitro_user_options_N.txt'
%
%   IPOPT options
%       60 - IPOPT_OPT, 0           See IPOPT_OPTIONS for details.
%
%   MINOPF options
%       61 - MNS_FEASTOL, 0 (1e-3)  primal feasibility tolerance,
%                                   set to value of OPF_VIOLATION by default
%       62 - MNS_ROWTOL, 0  (1e-3)  row tolerance
%                                   set to value of OPF_VIOLATION by default
%       63 - MNS_XTOL, 0    (1e-3)  x tolerance
%                                   set to value of CONSTR_TOL_X by default
%       64 - MNS_MAJDAMP, 0 (0.5)   major damping parameter
%       65 - MNS_MINDAMP, 0 (2.0)   minor damping parameter
%       66 - MNS_PENALTY_PARM, 0 (1.0)  penalty parameter
%       67 - MNS_MAJOR_IT, 0 (200)  major iterations
%       68 - MNS_MINOR_IT, 0 (2500) minor iterations
%       69 - MNS_MAX_IT, 0 (2500)   iterations limit
%       70 - MNS_VERBOSITY, -1
%           [  -1 - controlled by VERBOSE option                            ]
%           [   0 - print nothing                                           ]
%           [   1 - print only termination status message                   ]
%           [   2 - print termination status and screen progress            ]
%           [   3 - print screen progress, report file (usually fort.9)     ]
%       71 - MNS_CORE, 0 (1200 * nb + 2 * (nb + ng)^2) memory allocation
%       72 - MNS_SUPBASIC_LIM, 0 (2*nb + 2*ng) superbasics limit
%       73 - MNS_MULT_PRICE, 0 (30) multiple price
%
%   MIPS (including MIPS-sc), PDIPM, SC-PDIPM, and TRALM options
%       81 - PDIPM_FEASTOL, 0       feasibility (equality) tolerance
%                                   for MIPS, PDIPM and SC-PDIPM, set
%                                   to value of OPF_VIOLATION by default
%       82 - PDIPM_GRADTOL, 1e-6    gradient tolerance for MIPS, PDIPM
%                                   and SC-PDIPM
%       83 - PDIPM_COMPTOL, 1e-6    complementary condition (inequality)
%                                   tolerance for MIPS, PDIPM and SC-PDIPM
%       84 - PDIPM_COSTTOL, 1e-6    optimality tolerance for MIPS, PDIPM
%                                   and SC-PDIPM
%       85 - PDIPM_MAX_IT,  150     maximum number of iterations for MIPS,
%                                   PDIPM and SC-PDIPM
%       86 - SCPDIPM_RED_IT, 20     maximum number of MIPS-sc or SC-PDIPM
%                                   reductions per iteration
%       87 - TRALM_FEASTOL, 0       feasibility tolerance for TRALM
%                                   set to value of OPF_VIOLATION by default
%       88 - TRALM_PRIMETOL, 5e-4   primal variable tolerance for TRALM
%       89 - TRALM_DUALTOL, 5e-4    dual variable tolerance for TRALM
%       90 - TRALM_COSTTOL, 1e-5    optimality tolerance for TRALM
%       91 - TRALM_MAJOR_IT, 40     maximum number of major iterations
%       92 - TRALM_MINOR_IT, 100    maximum number of minor iterations
%       93 - SMOOTHING_RATIO, 0.04  piecewise linear curve smoothing ratio
%                                   used in SC-PDIPM and TRALM
%
%   CPLEX options
%       95 - CPLEX_LPMETHOD, 0      solution algorithm for continuous LPs
%           [   0 - automatic: let CPLEX choose                             ]
%           [   1 - primal simplex                                          ]
%           [   2 - dual simplex                                            ]
%           [   3 - network simplex                                         ]
%           [   4 - barrier                                                 ]
%           [   5 - sifting                                                 ]
%           [   6 - concurrent (dual, barrier, and primal)                  ]
%       96 - CPLEX_QPMETHOD, 0      solution algorithm for continuous QPs
%           [   0 - automatic: let CPLEX choose                             ]
%           [   1 - primal simplex optimizer                                ]
%           [   2 - dual simplex optimizer                                  ]
%           [   3 - network optimizer                                       ]
%           [   4 - barrier optimizer                                       ]
%       97 - CPLEX_OPT, 0           See CPLEX_OPTIONS for details
%
%   MOSEK options
%       111 - MOSEK_LP_ALG, 0       solution algorithm for continuous LPs
%                                       (MSK_IPAR_OPTIMIZER)
%           [   0 - automatic: let MOSEK choose                             ]
%           [   1 - interior point                                          ]
%           [   4 - primal simplex                                          ]
%           [   5 - dual simplex                                            ]
%           [   6 - primal dual simplex                                     ]
%           [   7 - automatic simplex (MOSEK chooses which simplex method)  ]
%           [   10 - concurrent                                             ]
%       112 - MOSEK_MAX_IT, 0 (400)     interior point max iterations
%                                           (MSK_IPAR_INTPNT_MAX_ITERATIONS)
%       113 - MOSEK_GAP_TOL, 0 (1e-8)   interior point relative gap tolerance
%                                           (MSK_DPAR_INTPNT_TOL_REL_GAP)
%       114 - MOSEK_MAX_TIME, 0 (-1)    maximum time allowed for solver
%                                           (MSK_DPAR_OPTIMIZER_MAX_TIME)
%       115 - MOSEK_NUM_THREADS, 0 (1)  maximum number of threads to use
%                                           (MSK_IPAR_INTPNT_NUM_THREADS)
%       116 - MOSEK_OPT, 0              See MOSEK_OPTIONS for details
%
%   Gurobi options
%       121 - GRB_METHOD, -1         solution algorithm (Method)
%           [  -1 - automatic, let Gurobi decide                            ]
%           [   0 - primal simplex                                          ]
%           [   1 - dual simplex                                            ]
%           [   2 - barrier                                                 ]
%           [   3 - concurrent (LP only)                                    ]
%           [   4 - deterministic concurrent (LP only)                      ]
%       122 - GRB_TIMELIMIT, Inf    maximum time allowed for solver (TimeLimit)
%       123 - GRB_THREADS, 0 (auto) maximum number of threads to use (Threads)
%       124 - GRB_OPT, 0            See GUROBI_OPTIONS for details

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%-----  set up default option values  -----
i = 1;
if rem(nargin, 2)       %% odd number of arguments
    options = varargin{1};  %% base options vector passed in
    i = 2;                  %% start processing parameters with 2nd one
else                    %% even number of parameters
    options = [             %% use defaults for base options vector
    
        %% power flow options
        1;      %% 1  - PF_ALG
        1e-8;   %% 2  - PF_TOL
        10;     %% 3  - PF_MAX_IT
        30;     %% 4  - PF_MAX_IT_FD
        1000;   %% 5  - PF_MAX_IT_GS
        0;      %% 6  - ENFORCE_Q_LIMS
        0;      %% 7  - RESERVED7
        0;      %% 8  - RESERVED8
        0;      %% 9  - RESERVED9
        0;      %% 10 - PF_DC
        
        %% OPF options
        0;      %% 11 - OPF_ALG
        0;      %% 12 - RESERVED12 (was OPF_ALG_POLY = 100)
        0;      %% 13 - RESERVED13 (was OPF_ALG_PWL = 200)
        0;      %% 14 - RESERVED14 (was OPF_POLY2PWL_PTS = 10)
        0;      %% 15 - OPF_NEQ (removed)
        5e-6;   %% 16 - OPF_VIOLATION
        1e-4;   %% 17 - CONSTR_TOL_X
        1e-4;   %% 18 - CONSTR_TOL_F
        0;      %% 19 - CONSTR_MAX_IT
        3e-3;   %% 20 - LPC_TOL_GRAD (removed)
        1e-4;   %% 21 - LPC_TOL_X (removed)
        400;    %% 22 - LPC_MAX_IT (removed)
        5;      %% 23 - LPC_MAX_RESTART (removed)
        0;      %% 24 - OPF_FLOW_LIM
        0;      %% 25 - OPF_IGNORE_ANG_LIM
        0;      %% 26 - OPF_ALG_DC
        0;      %% 27 - RESERVED27
        0;      %% 28 - RESERVED28
        0;      %% 29 - RESERVED29
        0;      %% 30 - RESERVED30
        
        %% output options
        1;      %% 31 - VERBOSE
        -1;     %% 32 - OUT_ALL
        1;      %% 33 - OUT_SYS_SUM
        0;      %% 34 - OUT_AREA_SUM
        1;      %% 35 - OUT_BUS
        1;      %% 36 - OUT_BRANCH
        0;      %% 37 - OUT_GEN
        -1;     %% 38 - OUT_ALL_LIM
        1;      %% 39 - OUT_V_LIM
        1;      %% 40 - OUT_LINE_LIM
        1;      %% 41 - OUT_PG_LIM
        1;      %% 42 - OUT_QG_LIM
        0;      %% 43 - RESERVED43 (was OUT_RAW)
        0;      %% 44 - OUT_FORCE
        0;      %% 45 - RESERVED45
        0;      %% 46 - RESERVED46
        0;      %% 47 - RESERVED47
        0;      %% 48 - RESERVED48
        0;      %% 49 - RESERVED49
        0;      %% 50 - RESERVED50
        
        %% other options
        1;      %% 51 - SPARSE_QP (removed)
        0;      %% 52 - RETURN_RAW_DER
        0;      %% 53 - RESERVED53
        0;      %% 54 - RESERVED54
        4;      %% 55 - FMC_ALG
        0;      %% 56 - RESERVED56
        0;      %% 57 - RESERVED57
        0;      %% 58 - KNITRO_OPT
        0;      %% 59 - RESERVED59
        0;      %% 60 - IPOPT_OPT
        
        %% MINOPF options
        0;      %% 61 - MNS_FEASTOL
        0;      %% 62 - MNS_ROWTOL
        0;      %% 63 - MNS_XTOL
        0;      %% 64 - MNS_MAJDAMP
        0;      %% 65 - MNS_MINDAMP
        0;      %% 66 - MNS_PENALTY_PARM
        0;      %% 67 - MNS_MAJOR_IT
        0;      %% 68 - MNS_MINOR_IT
        0;      %% 69 - MNS_MAX_IT
        -1;     %% 70 - MNS_VERBOSITY
        0;      %% 71 - MNS_CORE
        0;      %% 72 - MNS_SUPBASIC_LIM
        0;      %% 73 - MNS_MULT_PRICE
        0;      %% 74 - RESERVED74
        0;      %% 75 - RESERVED75
        0;      %% 76 - RESERVED76
        0;      %% 77 - RESERVED77
        0;      %% 78 - RESERVED78
        0;      %% 79 - RESERVED79
        0;      %% 80 - FORCE_PC_EQ_P0, for c3sopf
        
        %% MIPS, PDIPM, SC-PDIPM, and TRALM options
        0;      %% 81 - PDIPM_FEASTOL
        1e-6;   %% 82 - PDIPM_GRADTOL
        1e-6;   %% 83 - PDIPM_COMPTOL
        1e-6;   %% 84 - PDIPM_COSTTOL
        150;    %% 85 - PDIPM_MAX_IT
        20;     %% 86 - SCPDIPM_RED_IT
        0;      %% 87 - TRALM_FEASTOL
        5e-4;   %% 88 - TRALM_PRIMETOL
        5e-4;   %% 89 - TRALM_DUALTOL
        1e-5;   %% 90 - TRALM_COSTTOL
        40;     %% 91 - TRALM_MAJOR_IT
        100;    %% 92 - TRALM_MINOR_IT
        0.04;   %% 93 - SMOOTHING_RATIO
        0;      %% 94 - RESERVED94
        
        %% CPLEX options
        0;      %% 95 - CPLEX_LPMETHOD
        0;      %% 96 - CPLEX_QPMETHOD
        0;      %% 97 - CPLEX_OPT
        0;      %% 98 - RESERVED98
        0;      %% 99 - RESERVED99
        0;      %% 100 - RESERVED100
        0;      %% 101 - RESERVED101
        0;      %% 102 - RESERVED102
        0;      %% 103 - RESERVED103
        0;      %% 104 - RESERVED104
        0;      %% 105 - RESERVED105
        0;      %% 106 - RESERVED106
        0;      %% 107 - RESERVED107
        0;      %% 108 - RESERVED108
        0;      %% 109 - RESERVED109
        0;      %% 110 - RESERVED110

        %% MOSEK options
        0;      %% 111 - MOSEK_LP_ALG
        0;      %% 112 - MOSEK_MAX_IT
        0;      %% 113 - MOSEK_GAP_TOL
        0;      %% 114 - MOSEK_MAX_TIME
        0;      %% 115 - MOSEK_NUM_THREADS
        0;      %% 116 - MOSEK_OPT
        0;      %% 117 - RESERVED117
        0;      %% 118 - RESERVED118
        0;      %% 119 - RESERVED119
        0;      %% 120 - RESERVED120

        %% Gurobi options
        -1;     %% 121 - GRB_METHOD
        Inf;    %% 122 - GRB_TIMELIMIT
        0;      %% 123 - GRB_THREADS
        0;      %% 124 - GRB_OPT
    ];
end

%%-----  set up option names  -----
%% power flow options
names = char(   'PF_ALG', ...               %% 1
                'PF_TOL', ...               %% 2
                'PF_MAX_IT', ...            %% 3
                'PF_MAX_IT_FD', ...         %% 4
                'PF_MAX_IT_GS', ...         %% 5
                'ENFORCE_Q_LIMS', ...       %% 6
                'RESERVED7', ...            %% 7
                'RESERVED8', ...            %% 8
                'RESERVED9', ...            %% 9
                'PF_DC' );                  %% 10

%% OPF options
names = char(   names, ...
                'OPF_ALG', ...              %% 11
                'RESERVED12', ...           %% 12   (was OPF_ALG_POLY)
                'RESERVED13', ...           %% 13   (was OPF_ALG_PWL)
                'RESERVED14', ...           %% 14   (was OPF_POLY2PWL_PTS)
                'OPF_NEQ', ...              %% 15   (removed)
                'OPF_VIOLATION', ...        %% 16
                'CONSTR_TOL_X', ...         %% 17
                'CONSTR_TOL_F', ...         %% 18
                'CONSTR_MAX_IT', ...        %% 19
                'LPC_TOL_GRAD'  );          %% 20   (removed)
names = char(   names, ...
                'LPC_TOL_X', ...            %% 21   (removed)
                'LPC_MAX_IT', ...           %% 22   (removed)
                'LPC_MAX_RESTART', ...      %% 23   (removed)
                'OPF_FLOW_LIM', ...         %% 24
                'OPF_IGNORE_ANG_LIM', ...   %% 25
                'OPF_ALG_DC', ...           %% 26
                'RESERVED27', ...           %% 27
                'RESERVED28', ...           %% 28
                'RESERVED29', ...           %% 29
                'RESERVED30'    );          %% 30

%% output options
names = char(   names, ...
                'VERBOSE', ...              %% 31
                'OUT_ALL', ...              %% 32
                'OUT_SYS_SUM', ...          %% 33
                'OUT_AREA_SUM', ...         %% 34
                'OUT_BUS', ...              %% 35
                'OUT_BRANCH', ...           %% 36
                'OUT_GEN', ...              %% 37
                'OUT_ALL_LIM', ...          %% 38
                'OUT_V_LIM', ...            %% 39
                'OUT_LINE_LIM'  );          %% 40
names = char(   names, ...
                'OUT_PG_LIM', ...           %% 41
                'OUT_QG_LIM', ...           %% 42
                'RESERVED43', ...           %% 43 (was OUT_RAW)
                'OUT_FORCE', ...            %% 44
                'RESERVED45', ...           %% 45
                'RESERVED46', ...           %% 46
                'RESERVED47', ...           %% 47
                'RESERVED48', ...           %% 48
                'RESERVED49', ...           %% 49
                'RESERVED50'    );          %% 50
%% other options
names = char(   names, ...
                'SPARSE_QP', ...            %% 51   (removed)
                'RETURN_RAW_DER', ...       %% 52
                'RESERVED53', ...           %% 53
                'RESERVED54', ...           %% 54
                'FMC_ALG', ...              %% 55
                'RESERVED56', ...           %% 56
                'RESERVED57', ...           %% 57
                'KNITRO_OPT', ...           %% 58
                'RESERVED59', ...           %% 59
                'IPOPT_OPT'     );          %% 60
%% MINOS options
names = char(   names, ...
                'MNS_FEASTOL', ...          %% 61
                'MNS_ROWTOL', ...           %% 62
                'MNS_XTOL', ...             %% 63
                'MNS_MAJDAMP', ...          %% 64
                'MNS_MINDAMP', ...          %% 65
                'MNS_PENALTY_PARM', ...     %% 66
                'MNS_MAJOR_IT', ...         %% 67
                'MNS_MINOR_IT', ...         %% 68
                'MNS_MAX_IT', ...           %% 69
                'MNS_VERBOSITY' );          %% 70
%% other flags
names = char(   names, ...
                'MNS_CORE', ...             %% 71
                'MNS_SUPBASIC_LIM', ...     %% 72
                'MNS_MULT_PRICE', ...       %% 73
                'RESERVED74', ...           %% 74
                'RESERVED75', ...           %% 75
                'RESERVED76', ...           %% 76
                'RESERVED77', ...           %% 77
                'RESERVED78', ...           %% 78
                'RESERVED79', ...           %% 79
                'FORCE_PC_EQ_P0'    );      %% 80

%% MIPS, PDIPM, SC-PDIPM, and TRALM options
names = char(   names, ...
                'PDIPM_FEASTOL', ...        %% 81
                'PDIPM_GRADTOL', ...        %% 82
                'PDIPM_COMPTOL', ...        %% 83
                'PDIPM_COSTTOL', ...        %% 84
                'PDIPM_MAX_IT', ...         %% 85
                'SCPDIPM_RED_IT', ...       %% 86
                'TRALM_FEASTOL', ...        %% 87
                'TRALM_PRIMETOL', ...       %% 88
                'TRALM_DUALTOL', ...        %% 89
                'TRALM_COSTTOL', ...        %% 90
                'TRALM_MAJOR_IT', ...       %% 91
                'TRALM_MINOR_IT', ...       %% 92
                'SMOOTHING_RATIO'   );      %% 93

%% CPLEX options
names = char(   names, ...
                'RESERVED94', ...           %% 94
                'CPLEX_LPMETHOD', ...       %% 95
                'CPLEX_QPMETHOD', ...       %% 96
                'CPLEX_OPT', ...            %% 97
                'RESERVED98', ...           %% 98
                'RESERVED99', ...           %% 99
                'RESERVED100', ...          %% 100
                'RESERVED101', ...          %% 101
                'RESERVED102', ...          %% 102
                'RESERVED103', ...          %% 103
                'RESERVED104', ...          %% 104
                'RESERVED105', ...          %% 105
                'RESERVED106', ...          %% 106
                'RESERVED107', ...          %% 107
                'RESERVED108', ...          %% 108
                'RESERVED109', ...          %% 109
                'RESERVED110'   );          %% 110

%% MOSEK options
names = char(   names, ...
                'MOSEK_LP_ALG', ...         %% 111
                'MOSEK_MAX_IT', ...         %% 112
                'MOSEK_GAP_TOL', ...        %% 113
                'MOSEK_MAX_TIME', ...       %% 114
                'MOSEK_NUM_THREADS', ...    %% 115
                'MOSEK_OPT', ...            %% 116
                'RESERVED117', ...          %% 117
                'RESERVED118', ...          %% 118
                'RESERVED119', ...          %% 119
                'RESERVED120'   );          %% 120

%% Gurobi options
names = char(   names, ...
                'GRB_METHOD', ...           %% 121
                'GRB_TIMELIMIT', ...        %% 122
                'GRB_THREADS', ...          %% 123
                'GRB_OPT'   );              %% 124

%%-----  process parameters  -----
while i <= nargin
    %% get parameter name and value
    pname = varargin{i};
    pval  = varargin{i+1};
    
    %% get parameter index
    namestr = names';
    namestr = namestr(:)';
    namelen = size(names, 2);
    pidx = ceil(findstr([pname blanks(namelen-length(pname))], namestr) / namelen);
    if isempty(pidx)
        error('"%s" is not a valid named option', pname);
    end
    % fprintf('''%s'' (%d) = %d\n', pname, pidx, pval);

    %% update option
    options(pidx) = pval;

    i = i + 2;                              %% go to next parameter
end
