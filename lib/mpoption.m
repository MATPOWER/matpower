function opt = mpoption(varargin)
%MPOPTION  Used to set and retrieve a MATPOWER options struct.
%
%   OPT = MPOPTION
%       Returns the default options struct.
%
%   OPT = MPOPTION(OVERRIDES)
%       Returns the default options struct, with some fields overridden
%       by values from OVERRIDES, which can be a struct or the name of
%       a function that returns a struct.
%
%   OPT = MPOPTION(NAME1, VALUE1, NAME2, VALUE2, ...)
%       Same as previous, except override options are specified by NAME,
%       VALUE pairs. This can be used to set any part of the options
%       struct. The names can be individual fields or multi-level field
%       names with embedded periods. The values can be scalars or structs.
%
%       For backward compatibility, the NAMES and VALUES may correspond
%       to old-style MATPOWER option names (elements in the old-style
%       options vector) as well.
%
%   OPT = MPOPTION(OPT0)
%       Converts an old-style options vector OPT0 into the corresponding
%       options struct. If OPT0 is an options struct it does nothing.
%
%   OPT = MPOPTION(OPT0, OVERRIDES)
%       Applies overrides to an existing set of options, OPT0, which
%       can be an old-style options vector or an options struct.
%
%   OPT = MPOPTION(OPT0, NAME1, VALUE1, NAME2, VALUE2, ...)
%       Same as above except it uses the old-style options vector OPT0
%       as a base instead of the old default options vector.
%
%   OPT_VECTOR = MPOPTION(OPT, [])
%       Creates and returns an old-style options vector from an
%       options struct OPT.
%
%   Note: The use of old-style MATPOWER options vectors and their
%         names and values has been deprecated and will be removed
%         in a future version of MATPOWER. Until then, all uppercase
%         option names are not permitted for new top-level options.
%
%   Examples:
%       mpopt = mpoption('pf.alg', 'FDXB', 'pf.tol', 1e-4);
%       mpopt = mpoption(mpopt, 'opf.dc.solver', 'CPLEX', 'verbose', 2);
%
%The currently defined options are as follows:
%
%   name                    default     description [options]
%----------------------    ---------   ----------------------------------
%Model options:
%   model                   'AC'        AC vs. DC power flow model
%       [ 'AC' - use nonlinear AC model & corresponding algorithms/options  ]
%       [ 'DC' - use linear DC model & corresponding algorithms/options     ]
%
%Power Flow options:
%   pf.alg                  'NR'        AC power flow algorithm
%       [ 'NR'   - Newton's method                                          ]
%       [ 'FDXB' - Fast-Decoupled (XB version)                              ]
%       [ 'FDBX' - Fast-Decoupled (BX version)                              ]
%       [ 'GS'   - Gauss-Seidel                                             ]
%       [ 'PQSUM'- Power Summation method (radial networks only)            ]
%       [ 'ISUM' - Current Summation method (radial networks only)          ]
%       [ 'YSUM' - Admittance Summation method (radial networks only)       ]
%   pf.tol                  1e-8        termination tolerance on per unit
%                                       P & Q mismatch
%   pf.nr.max_it            10          maximum number of iterations for
%                                       Newton's method
%   pf.nr.lin_solver        ''          linear solver passed to MPLINSOLVE to
%                                       solve Newton update step
%       [ ''      - default to '\' for small systems, 'LU3' for larger ones ]
%       [ '\'     - built-in backslash operator                             ]
%       [ 'LU'    - explicit default LU decomposition and back substitution ]
%       [ 'LU3'   - 3 output arg form of LU, Gilbert-Peierls algorithm with ]
%       [           approximate minimum degree (AMD) reordering             ]
%       [ 'LU4'   - 4 output arg form of LU, UMFPACK solver (same as 'LU')  ]
%       [ 'LU5'   - 5 output arg form of LU, UMFPACK solver, w/row scaling  ]
%       [       (see MPLINSOLVE for complete list of all options)           ]
%   pf.fd.max_it            30          maximum number of iterations for
%                                       fast decoupled method
%   pf.gs.max_it            1000        maximum number of iterations for
%                                       Gauss-Seidel method
%   pf.radial.max_it        20          maximum number of iterations for
%                                       radial power flow methods
%   pf.radial.vcorr         0           perform voltage correction procedure
%                                       in distribution power flow
%       [  0 - do NOT perform voltage correction                            ]
%       [  1 - perform voltage correction                                   ]
%   pf.enforce_q_lims       0           enforce gen reactive power limits at
%                                       expense of |V|
%       [  0 - do NOT enforce limits                                        ]
%       [  1 - enforce limits, simultaneous bus type conversion             ]
%       [  2 - enforce limits, one-at-a-time bus type conversion            ]
%
%Continuation Power Flow options:
%   cpf.parameterization    3           choice of parameterization
%       [  1 - natural                                                      ]
%       [  2 - arc length                                                   ]
%       [  3 - pseudo arc length                                            ]
%   cpf.stop_at             'NOSE'      determins stopping criterion
%       [ 'NOSE'     - stop when nose point is reached                      ]
%       [ 'FULL'     - trace full nose curve                                ]
%       [ <lam_stop> - stop upon reaching specified target lambda value     ]
%   cpf.enforce_p_lims      0           enforce gen active power limits
%       [  0 - do NOT enforce limits                                        ]
%       [  1 - enforce limits, simultaneous bus type conversion             ]
%   cpf.enforce_q_lims      0           enforce gen reactive power limits at
%                                       expense of |V|
%       [  0 - do NOT enforce limits                                        ]
%       [  1 - enforce limits, simultaneous bus type conversion             ]
%   cpf.enforce_v_lims      0           enforce bus voltage magnitude limits
%       [  0 - do NOT enforce limits                                        ]
%       [  1 - enforce limits, termination on detection                     ]
%   cpf.enforce_flow_lims   0           enforce branch flow MVA limits
%       [  0 - do NOT enforce limits                                        ]
%       [  1 - enforce limits, termination on detection                     ]
%   cpf.step                0.05        continuation power flow step size
%   cpf.adapt_step          0           toggle adaptive step size feature
%       [  0 - adaptive step size disabled                                  ]
%       [  1 - adaptive step size enabled                                   ]
%   cpf.step_min            1e-4        minimum allowed step size
%   cpf.step_max            0.2         maximum allowed step size
%   cpf.adapt_step_damping  0.7         damping factor for adaptive step
%                                       sizing
%   cpf.adapt_step_tol      1e-3        tolerance for adaptive step sizing
%   cpf.target_lam_tol      1e-5        tolerance for target lambda detection
%   cpf.nose_tol            1e-5        tolerance for nose point detection (pu)
%   cpf.p_lims_tol          0.01        tolerance for generator active
%                                       power limit enforcement (MW)
%   cpf.q_lims_tol          0.01        tolerance for generator reactive
%                                       power limit enforcement (MVAR)
%   cpf.v_lims_tol          1e-4        tolerance for bus voltage
%                                       magnitude enforcement (p.u)
%   cpf.flow_lims_tol       0.01        tolerance for line MVA flow
%                                       enforcement (MVA)
%   cpf.plot.level          0           control plotting of noze curve
%       [  0 - do not plot nose curve                                       ]
%       [  1 - plot when completed                                          ]
%       [  2 - plot incrementally at each iteration                         ]
%       [  3 - same as 2, with 'pause' at each iteration                    ]
%   cpf.plot.bus            <empty>     index of bus whose voltage is to be
%                                       plotted
%   cpf.user_callback       <empty>     string containing the name of a user
%                                       callback function, or struct with
%                                       function name, and optional priority
%                                       and/or args, or cell array of such
%                                       strings and/or structs, see
%                                       'help cpf_default_callback' for details
%
%Optimal Power Flow options:
%   name                    default     description [options]
%----------------------    ---------   ----------------------------------
%   opf.ac.solver           'DEFAULT'   AC optimal power flow solver
%       [ 'DEFAULT' - choose default solver, i.e. 'MIPS'                    ]
%       [ 'MIPS'    - MIPS, MATPOWER Interior Point Solver, primal/dual     ]
%       [             interior point method (pure MATLAB/Octave)            ]
%       [ 'FMINCON' - MATLAB Optimization Toolbox, FMINCON                  ]
%       [ 'IPOPT'   - IPOPT, requires MEX interface to IPOPT solver         ]
%       [             available from:                                       ]
%       [                 http://www.coin-or.org/projects/Ipopt.xml         ]
%       [ 'KNITRO'  - KNITRO, requires MATLAB Optimization Toolbox and      ]
%       [             KNITRO libraries available from: http://www.ziena.com/]
%       [ 'MINOPF'  - MINOPF, MINOS-based solver, requires optional         ]
%       [             MEX-based MINOPF package, available from:             ]
%       [                   http://www.pserc.cornell.edu/minopf/            ]
%       [ 'PDIPM'   - PDIPM, primal/dual interior point method, requires    ]
%       [             optional MEX-based TSPOPF package, available from:    ]
%       [                   http://www.pserc.cornell.edu/tspopf/            ]
%       [ 'SDPOPF'  - SDPOPF, solver based on semidefinite relaxation of    ]
%       [             OPF problem, requires optional packages:              ]
%       [               SDP_PF, available in extras/sdp_pf                  ]
%       [               YALMIP, available from:                             ]
%       [                   http://users.isy.liu.se/johanl/yalmip/          ]
%       [               SDP solver such as SeDuMi, available from:          ]
%       [                   http://sedumi.ie.lehigh.edu/                    ]
%       [ 'TRALM'   - TRALM, trust region based augmented Langrangian       ]
%       [             method, requires TSPOPF (see 'PDIPM')                 ]
%   opf.dc.solver           'DEFAULT'   DC optimal power flow solver
%       [ 'DEFAULT' - choose solver based on availability in the following  ]
%       [             order: 'GUROBI', 'CPLEX', 'MOSEK', 'OT',              ]
%       [             'GLPK' (linear costs only), 'BPMPD', 'MIPS'           ]
%       [ 'MIPS'    - MIPS, MATPOWER Interior Point Solver, primal/dual     ]
%       [             interior point method (pure MATLAB/Octave)            ]
%       [ 'BPMPD'   - BPMPD, requires optional MEX-based BPMPD_MEX package  ]
%       [             available from: http://www.pserc.cornell.edu/bpmpd/   ]
%       [ 'CLP'     - CLP, requires interface to COIN-OP LP solver          ]
%       [             available from:http://www.coin-or.org/projects/Clp.xml]
%       [ 'CPLEX'   - CPLEX, requires CPLEX solver available from:          ]
%       [             http://www.ibm.com/software/integration/ ...          ]
%       [                                 ... optimization/cplex-optimizer/ ]
%       [ 'GLPK'    - GLPK, requires interface to GLPK solver               ]
%       [             available from: http://www.gnu.org/software/glpk/     ]
%       [             (GLPK does not work with quadratic cost functions)    ]
%       [ 'GUROBI'  - GUROBI, requires Gurobi optimizer (v. 5+)             ]
%       [             available from: http://www.gurobi.com/                ]
%       [ 'IPOPT'   - IPOPT, requires MEX interface to IPOPT solver         ]
%       [             available from:                                       ]
%       [                 http://www.coin-or.org/projects/Ipopt.xml         ]
%       [ 'MOSEK'   - MOSEK, requires MATLAB interface to MOSEK solver      ]
%       [             available from: http://www.mosek.com/                 ]
%       [ 'OT'      - MATLAB Optimization Toolbox, QUADPROG, LINPROG        ]
%   opf.current_balance     0           type of nodal balance equation
%       [ 0 - use complex power balance equations                           ]
%       [ 1 - use complex current balance equations                         ]
%   opf.v_cartesian         0           voltage representation
%       [ 0 - bus voltage variables represented in polar coordinates        ]
%       [ 1 - bus voltage variables represented in Cartesian coordinates    ]
%   opf.violation           5e-6        constraint violation tolerance
%   opf.use_vg              0           respect gen voltage setpt     [ 0-1 ]
%       [ 0 - use specified bus Vmin & Vmax, and ignore gen Vg              ]
%       [ 1 - replace specified bus Vmin & Vmax by corresponding gen Vg     ]
%       [ between 0 and 1 - use a weighted average of the 2 options         ]
%   opf.flow_lim            'S'         quantity limited by branch flow
%                                       constraints
%       [ 'S' - apparent power flow (limit in MVA)                          ]
%       [ 'P' - active power flow, implemented using P   (limit in MW)      ]
%       [ '2' - active power flow, implemented using P^2 (limit in MW)      ]
%       [ 'I' - current magnitude (limit in MVA at 1 p.u. voltage)          ]
%   opf.ignore_angle_lim    0           angle diff limits for branches
%       [ 0 - include angle difference limits, if specified                 ]
%       [ 1 - ignore angle difference limits even if specified              ]
%   opf.softlims.default    1           behavior of OPF soft limits for
%                                       which parameters are not explicitly
%                                       provided
%       [ 0 - do not include softlims if not explicitly specified           ]
%       [ 1 - include softlims w/default values if not explicitly specified ]
%   opf.init_from_mpc       -1          (DEPRECATED: use opf.start instead)
%                                       specify whether to use current state
%                                       in MATPOWER case to initialize OPF
%                                       (currently only supported by fmincon,
%                                        Ipopt, Knitro and MIPS solvers,
%                                        others always use mpc)
%       [  -1 - MATPOWER decides, based on solver/algorithm                 ]
%       [   0 - ignore current state in MATPOWER case (only applies to      ]
%       [       fmincon, Ipopt, Knitro and MIPS, which use an interior pt   ]
%       [       estimate; others use current state as with opf.start = 2)   ]
%       [   1 - use current state in MATPOWER case                          ]
%   opf.start               0           strategy for initializing OPF start pt
%       [   0 - default, MATPOWER decides based on solver                   ]
%       [       (currently identical to 1)                                  ]
%       [   1 - ignore current state in MATPOWER case (only applies to      ]
%       [       fmincon, Ipopt, Knitro and MIPS, which use an interior pt   ]
%       [       estimate; others use current state as with opf.start = 2)   ]
%       [   2 - use current state in MATPOWER case                          ]
%       [   3 - solve power flow and use resulting state                    ]
%   opf.return_raw_der      0           for AC OPF, return constraint and
%                                       derivative info in results.raw
%                                       (in fields g, dg, df, d2f) [ 0 or 1 ]
%
%Output options:
%   name                    default     description [options]
%----------------------    ---------   ----------------------------------
%   verbose                 1           amount of progress info printed
%       [   0 - print no progress info                                      ]
%       [   1 - print a little progress info                                ]
%       [   2 - print a lot of progress info                                ]
%       [   3 - print all progress info                                     ]
%   out.all                 -1          controls pretty-printing of results
%       [  -1 - individual flags control what prints                        ]
%       [   0 - do not print anything (overrides individual flags, ignored  ]
%       [       for files specified as FNAME arg to runpf(), runopf(), etc.)]
%       [   1 - print everything (overrides individual flags)               ]
%   out.sys_sum             1           print system summary       [ 0 or 1 ]
%   out.area_sum            0           print area summaries       [ 0 or 1 ]
%   out.bus                 1           print bus detail           [ 0 or 1 ]
%   out.branch              1           print branch detail        [ 0 or 1 ]
%   out.gen                 0           print generator detail     [ 0 or 1 ]
%   out.lim.all             -1          controls constraint info output
%       [  -1 - individual flags control what constraint info prints        ]
%       [   0 - no constraint info (overrides individual flags)             ]
%       [   1 - binding constraint info (overrides individual flags)        ]
%       [   2 - all constraint info (overrides individual flags)            ]
%   out.lim.v               1           control voltage limit info
%       [   0 - do not print                                                ]
%       [   1 - print binding constraints only                              ]
%       [   2 - print all constraints                                       ]
%       [   (same options for OUT_LINE_LIM, OUT_PG_LIM, OUT_QG_LIM)         ]
%   out.lim.line            1           control line flow limit info
%   out.lim.pg              1           control gen active power limit info
%   out.lim.qg              1           control gen reactive pwr limit info
%   out.force               0           print results even if success
%                                       flag = 0                   [ 0 or 1 ]
%   out.suppress_detail     -1          suppress all output but system summary
%       [  -1 - suppress details for large systems (> 500 buses)            ]
%       [   0 - do not suppress any output specified by other flags         ]
%       [   1 - suppress all output except system summary section           ]
%       [       (overrides individual flags, but not out.all = 1)           ]
%
%Solver specific options:
%       name                    default     description [options]
%   -----------------------    ---------   ----------------------------------
%   MIPS:
%       mips.linsolver          ''          linear system solver
%           [   '' or '\'   build-in backslash \ operator (e.g. x = A \ b)  ]
%           [   'PARDISO'   PARDISO solver (if available)                   ]
%       mips.feastol            0           feasibility (equality) tolerance
%                                           (set to opf.violation by default)
%       mips.gradtol            1e-6        gradient tolerance
%       mips.comptol            1e-6        complementary condition
%                                           (inequality) tolerance
%       mips.costtol            1e-6        optimality tolerance
%       mips.max_it             150         maximum number of iterations
%       mips.step_control       0           enable step-size cntrl [ 0 or 1 ]
%       mips.sc.red_it          20          maximum number of reductions per
%                                           iteration with step control
%       mips.xi                 0.99995     constant used in alpha updates*
%       mips.sigma              0.1         centering parameter*
%       mips.z0                 1           used to initialize slack variables*
%       mips.alpha_min          1e-8        returns "Numerically Failed" if
%                                           either alpha parameter becomes
%                                           smaller than this value*
%       mips.rho_min            0.95        lower bound on rho_t*
%       mips.rho_max            1.05        upper bound on rho_t*
%       mips.mu_threshold       1e-5        KT multipliers smaller than this
%                                           value for non-binding constraints
%                                           are forced to zero
%       mips.max_stepsize       1e10        returns "Numerically Failed" if the
%                                           2-norm of the reduced Newton step
%                                           exceeds this value*
%           * See the corresponding Appendix in the manual for details.
%
%   CPLEX:
%       cplex.lpmethod          0           solution algorithm for LP problems
%           [   0 - automatic: let CPLEX choose                             ]
%           [   1 - primal simplex                                          ]
%           [   2 - dual simplex                                            ]
%           [   3 - network simplex                                         ]
%           [   4 - barrier                                                 ]
%           [   5 - sifting                                                 ]
%           [   6 - concurrent (dual, barrier, and primal)                  ]
%       cplex.qpmethod          0           solution algorithm for QP problems
%           [   0 - automatic: let CPLEX choose                             ]
%           [   1 - primal simplex optimizer                                ]
%           [   2 - dual simplex optimizer                                  ]
%           [   3 - network optimizer                                       ]
%           [   4 - barrier optimizer                                       ]
%       cplex.opts              <empty>     see CPLEX_OPTIONS for details
%       cplex.opt_fname         <empty>     see CPLEX_OPTIONS for details
%       cplex.opt               0           see CPLEX_OPTIONS for details
%
%   FMINCON:
%       fmincon.alg             4           algorithm used by fmincon() for OPF
%                                           for Opt Toolbox 4 and later
%            [  1 - active-set (not suitable for large problems)            ]
%            [  2 - interior-point, w/default 'bfgs' Hessian approx         ]
%            [  3 - interior-point, w/ 'lbfgs' Hessian approx               ]
%            [  4 - interior-point, w/exact user-supplied Hessian           ]
%            [  5 - interior-point, w/Hessian via finite differences        ]
%            [  6 - sqp (not suitable for large problems)                   ]
%       fmincon.tol_x           1e-4        termination tol on x
%       fmincon.tol_f           1e-4        termination tol on f
%       fmincon.max_it          0           maximum number of iterations
%                                                           [  0 => default ]
%
%   GUROBI:
%       gurobi.method           0           solution algorithm (Method)
%           [  -1 - automatic, let Gurobi decide                            ]
%           [   0 - primal simplex                                          ]
%           [   1 - dual simplex                                            ]
%           [   2 - barrier                                                 ]
%           [   3 - concurrent (LP only)                                    ]
%           [   4 - deterministic concurrent (LP only)                      ]
%       gurobi.timelimit        Inf         maximum time allowed (TimeLimit)
%       gurobi.threads          0           max number of threads (Threads)
%       gurobi.opts             <empty>     see GUROBI_OPTIONS for details
%       gurobi.opt_fname        <empty>     see GUROBI_OPTIONS for details
%       gurobi.opt              0           see GUROBI_OPTIONS for details
%
%   IPOPT:
%       ipopt.opts              <empty>     see IPOPT_OPTIONS for details
%       ipopt.opt_fname         <empty>     see IPOPT_OPTIONS for details
%       ipopt.opt               0           see IPOPT_OPTIONS for details
%
%   KNITRO:
%       knitro.tol_x            1e-4        termination tol on x
%       knitro.tol_f            1e-4        termination tol on f
%       knitro.maxit            0           maximum number of iterations
%                                                           [  0 => default ]
%       knitro.opt_fname        <empty>     name of user-supplied native
%                                           KNITRO options file that overrides
%                                           all other options
%       knitro.opt              0           if knitro.opt_fname is empty and
%                                           knitro.opt is a non-zero integer N
%                                           then knitro.opt_fname is auto-
%                                           generated as:
%                                           'knitro_user_options_N.txt'
%
%   LINPROG:
%       linprog                 <empty>     LINPROG options passed to
%                                           OPTIMOPTIONS or OPTIMSET.
%                                           see LINPROG in the Optimization
%                                           Toolbox for details
%
%   MINOPF:
%       minopf.feastol          0 (1e-3)    primal feasibility tolerance
%                                           (set to opf.violation by default)
%       minopf.rowtol           0 (1e-3)    row tolerance
%       minopf.xtol             0 (1e-4)    x tolerance
%       minopf.majdamp          0 (0.5)     major damping parameter
%       minopf.mindamp          0 (2.0)     minor damping parameter
%       minopf.penalty          0 (1.0)     penalty parameter
%       minopf.major_it         0 (200)     major iterations
%       minopf.minor_it         0 (2500)    minor iterations
%       minopf.max_it           0 (2500)    iterations limit
%       minopf.verbosity        -1          amount of progress info printed
%           [  -1 - controlled by 'verbose' option                          ]
%           [   0 - print nothing                                           ]
%           [   1 - print only termination status message                   ]
%           [   2 - print termination status and screen progress            ]
%           [   3 - print screen progress, report file (usually fort.9)     ]
%       minopf.core             0 (1200*nb + 2*(nb+ng)^2) memory allocation
%       minopf.supbasic_lim     0 (2*nb + 2*ng) superbasics limit
%       minopf.mult_price       0 (30)      multiple price
%
%   MOSEK:
%       mosek.lp_alg            0           solution algorithm
%                                               (MSK_IPAR_OPTIMIZER)
%           for MOSEK 8.x ...        (see MOSEK_SYMBCON for a "better way")
%           [   0 - automatic: let MOSEK choose                             ]
%           [   1 - dual simplex                                            ]
%           [   2 - automatic: let MOSEK choose                             ]
%           [   3 - automatic simplex (MOSEK chooses which simplex method)  ]
%           [   4 - interior point                                          ]
%           [   6 - primal simplex                                          ]
%       mosek.max_it            0 (400)     interior point max iterations
%                                               (MSK_IPAR_INTPNT_MAX_ITERATIONS)
%       mosek.gap_tol           0 (1e-8)    interior point relative gap tol
%                                               (MSK_DPAR_INTPNT_TOL_REL_GAP)
%       mosek.max_time          0 (-1)      maximum time allowed
%                                               (MSK_DPAR_OPTIMIZER_MAX_TIME)
%       mosek.num_threads       0 (1)       max number of threads
%                                               (MSK_IPAR_INTPNT_NUM_THREADS)
%       mosek.opts              <empty>     see MOSEK_OPTIONS for details
%       mosek.opt_fname         <empty>     see MOSEK_OPTIONS for details
%       mosek.opt               0           see MOSEK_OPTIONS for details
%
%   QUADPROG:
%       quadprog                <empty>     QUADPROG options passed to
%                                           OPTIMOPTIONS or OPTIMSET.
%                                           see QUADPROG in the Optimization
%                                           Toolbox for details
%
%   TSPOPF:
%       pdipm.feastol           0           feasibility (equality) tolerance
%                                           (set to opf.violation by default)
%       pdipm.gradtol           1e-6        gradient tolerance
%       pdipm.comptol           1e-6        complementary condition
%                                           (inequality) tolerance
%       pdipm.costtol           1e-6        optimality tolerance
%       pdipm.max_it            150         maximum number of iterations
%       pdipm.step_control      0           enable step-size cntrl [ 0 or 1 ]
%       pdipm.sc.red_it         20          maximum number of reductions per
%                                           iteration with step control
%       pdipm.sc.smooth_ratio   0.04        piecewise linear curve smoothing
%                                           ratio
%
%       tralm.feastol           0           feasibility tolerance
%                                           (set to opf.violation by default)
%       tralm.primaltol         5e-4        primal variable tolerance
%       tralm.dualtol           5e-4        dual variable tolerance
%       tralm.costtol           1e-5        optimality tolerance
%       tralm.major_it          40          maximum number of major iterations
%       tralm.minor_it          40          maximum number of minor iterations
%       tralm.smooth_ratio      0.04        piecewise linear curve smoothing
%                                           ratio
%
%Experimental Options:
%   exp.sys_wide_zip_loads.pw   <empty>     1 x 3 vector of active load fraction
%                                           to be modeled as constant power,
%                                           constant current and constant
%                                           impedance, respectively, where
%                                           <empty> means use [1 0 0]
%   exp.sys_wide_zip_loads.qw   <empty>     same for reactive power, where
%                                           <empty> means use same value as
%                                           for 'pw'

%   MATPOWER
%   Copyright (c) 2013-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% some constants
N = 124;                %% number of options in old-style vector (MATPOWER 4.1)
N40 = 116;              %% dimension of MATPOWER 4.0 options vector
N32 = 93;               %% dimension of MATPOWER 3.2 options vector
v = mpoption_version;   %% version number of MATPOWER options struct

%% initialize flags and arg counter
have_opt0 = 0;          %% existing options struct or vector provided?
have_old_style_ov = 0;  %% override options using old-style names?
return_old_style = 0;   %% return value as old-style vector?
k = 1;
if nargin > 0
    opt0 = varargin{k};
    if isstruct(opt0) && isfield(opt0, 'v')         %% options struct
        have_opt0 = 1;
        k = k + 1;
    elseif isnumeric(opt0) && size(opt0, 2) == 1    %% options vector
        nn = size(opt0, 1);
        if ismember(nn, [N N40 N32])                %% of valid size
            %% expand smaller option vectors (from before MATPOWER 4.1)
            if nn < N
                optv = mpoption_old();
                opt0(nn+1:N) = optv(nn+1:N);
            end
            have_opt0 = 1;
            k = k + 1;
        end
    end
end

%% create base options vector to which overrides are made
if have_opt0
    if isstruct(opt0)               %% it's already a valid options struct
        if DEBUG, fprintf('OPT0 is a valid options struct\n'); end
        if opt0.v < v
            %% convert older version to current version
            opt_d = mpoption_default();
            if opt0.v == 1          %% convert version 1 to 2
                if isfield(opt_d, 'linprog')
                    opt0.lingprog = opt_d.linprog;
                end
                if isfield(opt_d, 'quadprog')
                    opt0.quadprog = opt_d.quadprog;
                end
            end
            if opt0.v <= 2          %% convert version 2 to 3
                opt0.out.suppress_detail = opt_d.out.suppress_detail;
            end
            %if opt0.v <= 3          %% convert version 3 to 4
                %% new mips options were all optional, no conversion needed
            %end
            if opt0.v <= 4          %% convert version 4 to 5
                opt0.opf.init_from_mpc = opt_d.opf.init_from_mpc;
            end
            if opt0.v <= 5          %% convert version 5 to 6
                if isfield(opt_d, 'clp')
                    opt0.clp = opt_d.clp;
                end
            end
            if opt0.v <= 6          %% convert version 6 to 7
                if isfield(opt_d, 'intlinprog')
                    opt0.intlinprog = opt_d.intlinprog;
                end
            end
            if opt0.v <= 7          %% convert version 7 to 8
                opt0.mips.linsolver = opt_d.mips.linsolver;
            end
            if opt0.v <= 8          %% convert version 8 to 9
                opt0.exp.sys_wide_zip_loads = opt_d.exp.sys_wide_zip_loads;
            end
            if opt0.v <= 9          %% convert version 9 to 10
                opt0.most = opt_d.most;
            end
            if opt0.v <= 10         %% convert version 10 to 11
                opt0.cpf.enforce_p_lims     = opt_d.cpf.enforce_p_lims;
                opt0.cpf.enforce_q_lims     = opt_d.cpf.enforce_q_lims;
                opt0.cpf.adapt_step_damping = opt_d.cpf.adapt_step_damping;
                opt0.cpf.target_lam_tol     = opt_d.cpf.target_lam_tol;
                opt0.cpf.nose_tol           = opt_d.cpf.nose_tol;
                opt0.cpf.p_lims_tol         = opt_d.cpf.p_lims_tol;
                opt0.cpf.q_lims_tol         = opt_d.cpf.q_lims_tol;
                if (~isempty(opt0.cpf.user_callback_args) && ...
                        ~isstruct(opt0.cpf.user_callback_args)) || ...
                        (isstruct(opt0.cpf.user_callback_args) && ...
                         ~isempty(fields(opt0.cpf.user_callback_args)))
                    warning('The ''cpf.user_callback_args'' option has been removed. Please include the args in a struct in ''cpf.user_callback'' instead.')
                end
                opt0.cpf = rmfield(opt0.cpf, 'user_callback_args');
            end
            if opt0.v <= 11         %% convert version 11 to 12
                opt0.opf.use_vg = opt_d.opf.use_vg;
            end
            if opt0.v <= 12         %% convert version 12 to 13
                opt0.pf.radial.max_it   = opt_d.pf.radial.max_it;
                opt0.pf.radial.vcorr    = opt_d.pf.radial.vcorr;
            end
            if opt0.v <= 13         %% convert version 13 to 14
                opt0.pf.nr.lin_solver   = opt_d.pf.nr.lin_solver;
            end
            if opt0.v <= 14         %% convert version 14 to 15
                opt0.cpf.enforce_v_lims     = opt_d.cpf.enforce_v_lims;
                opt0.cpf.enforce_flow_lims  = opt_d.cpf.enforce_flow_lims;
                opt0.cpf.v_lims_tol         = opt_d.cpf.v_lims_tol;
                opt0.cpf.flow_lims_tol      = opt_d.cpf.flow_lims_tol;
            end
            if opt0.v <= 15          %% convert version 15 to 16
                opt0.opf.start  = opt_d.opf.start;
            end
            if opt0.v <= 16          %% convert version 16 to 17
                if isfield(opt_d, 'knitro')
                    opt0.knitro.maxit  = opt_d.knitro.maxit;
                end
            end
            if opt0.v <= 17         %% convert version 17 to 18
                opt0.opf.current_balance = opt_d.opf.current_balance;
                opt0.opf.v_cartesian     = opt_d.opf.v_cartesian;
            end
            if opt0.v <= 18         %% convert version 18 to 19
                opt0.opf.softlims.default = opt_d.opf.softlims.default;
            end
            opt0.v = v;
        end
        opt = opt0;
    else                            %% convert from old-style options vector
        if DEBUG, fprintf('OPT0 is a old-style options vector\n'); end
        opt = mpoption_v2s(opt0);
    end
else                                %% use default options struct as base
    if DEBUG, fprintf('no OPT0, starting with default options struct\n'); end
    opt = mpoption_default();
end


%% do we have OVERRIDES or NAME/VALUE pairs
ov = [];
if nargin - k == 0          %% looking at last arg, must be OVERRIDES
    if isstruct(varargin{k})        %% OVERRIDES provided as struct
        if DEBUG, fprintf('OVERRIDES struct\n'); end
        ov = varargin{k};
    elseif ischar(varargin{k})      %% OVERRIDES provided as file/function name
        if DEBUG, fprintf('OVERRIDES file/function name\n'); end
        try
            ov = feval(varargin{k});
        catch
            error('mpoption: Unable to load MATPOWER options from ''%s''', varargin{k});
        end
        if ~isstruct(ov)
            error('mpoption: calling ''%s'' did not return a struct', varargin{k});
        end
    elseif isempty(varargin{k})
        return_old_style = 1;
    else
        error('mpoption: OVERRIDES must be a struct or the name of a function that returns a struct');
    end
elseif nargin - k > 0 && mod(nargin-k, 2)   %% even number of remaining args
    if DEBUG, fprintf('NAME/VALUE pairs override defaults\n'); end
    %% process NAME/VALUE pairs
    if strcmp(varargin{k}, upper(varargin{k}))  %% old-style, all UPPERCASE option pairs
        %% NOTE: new top-level option fields cannot be all uppercase
        if ~have_opt0
            opt_v = mpoption_old(varargin{:});  %% create modified vector ...
            opt = mpoption_v2s(opt_v);          %% ... then convert
        else
            have_old_style_ov = 1;
            %% convert pairs to struct
            while k < nargin
                name = varargin{k};
                val  = varargin{k+1};
                k = k + 2;
                ov.(name) = val;
            end
        end
    else                                        %% new option pairs
        %% convert pairs to struct
        while k < nargin
            name = varargin{k};
            val  = varargin{k+1};
            k = k + 2;
            c = regexp(name, '([^\.]*)', 'tokens');
            s = struct();
            for i = 1:length(c)
                s(i).type = '.';
                s(i).subs = c{i}{1};
            end
            ov = subsasgn(ov, s, val);
        end
    end
elseif nargin == 0 || nargin == 1
    if DEBUG, fprintf('no OVERRIDES, return default options struct or converted OPT0 vector\n'); end
else
    error('mpoption: invalid calling syntax, see ''help mpoption'' to double-check the valid options');
end

%% apply overrides
if ~isempty(ov)
    if have_old_style_ov
        opt = apply_old_mpoption_overrides(opt, ov);
    else
        persistent nsc_opt;     %% cache this to speed things up
        if ~isstruct(nsc_opt)
            vf = nested_struct_copy(mpoption_default(), mpoption_info_mips('V'));
            vf = nested_struct_copy(vf, mpoption_optional_fields());
            ex = struct(...
                'name', {}, ...
                'check', {}, ...
                'copy_mode', {} ...
            );
            %% add exceptions for optional packages
            opt_pkgs = mpoption_optional_pkgs();
            n = length(ex);
            for k = 1:length(opt_pkgs)
                fname = ['mpoption_info_' opt_pkgs{k}];
                if exist(fname, 'file') == 2
                    opt_ex = feval(fname, 'E');
                    nex = length(opt_ex);
                    if ~isempty(opt_ex)
                        for j = 1:nex
                            ex(n+j).name = opt_ex(j).name;
                        end
                        if isfield(opt_ex, 'check')
                            for j = 1:nex
                                ex(n+j).check = opt_ex(j).check;
                            end
                        end
                        if isfield(opt_ex, 'copy_mode')
                            for j = 1:nex
                                ex(n+j).copy_mode = opt_ex(j).copy_mode;
                            end
                        end
                        if isfield(opt_ex, 'valid_fields')
                            for j = 1:nex
                                ex(n+j).valid_fields = opt_ex(j).valid_fields;
                            end
                        end
                        n = n + nex;
                    end
                end
            end
            nsc_opt = struct('check', 1, 'valid_fields', vf, 'exceptions', ex);
        end
%         if have_fcn('catchme')
%             try
%                 opt = nested_struct_copy(opt, ov, nsc_opt);
%             catch me
%                 str = strrep(me.message, 'field', 'option');
%                 str = strrep(str, 'nested_struct_copy', 'mpoption');
%                 error(str);
%             end
%         else
            try
                opt = nested_struct_copy(opt, ov, nsc_opt);
            catch
                me = lasterr;
                str = strrep(me, 'field', 'option');
                str = strrep(str, 'nested_struct_copy', 'mpoption');
                error(str);
            end
%         end
    end
end
if return_old_style
    opt = mpoption_s2v(opt);
end


%%-------------------------------------------------------------------
function opt = apply_old_mpoption_overrides(opt0, ov)
%
%   OPT0 is assumed to already have all of the fields and sub-fields found
%   in the default options struct.

%% initialize output
opt = opt0;

errstr = 'mpoption: %g is not a valid value for the old-style ''%s'' option';
fields = fieldnames(ov);
for f = 1:length(fields)
    ff = fields{f};
    switch ff
        case 'PF_ALG'
            switch ov.(ff)
                case 1
                    opt.pf.alg = 'NR';      %% Newton's method
                case 2
                    opt.pf.alg = 'FDXB';    %% fast-decoupled (XB version)
                case 3
                    opt.pf.alg = 'FDBX';    %% fast-decoupled (BX version)
                case 4
                    opt.pf.alg = 'GS';      %% Gauss-Seidel
                otherwise
                    error(errstr, ov.(ff), ff);
            end
        case 'PF_TOL'
            opt.pf.tol = ov.(ff);
        case 'PF_MAX_IT'
            opt.pf.nr.max_it = ov.(ff);
        case 'PF_MAX_IT_FD'
            opt.pf.fd.max_it = ov.(ff);
        case 'PF_MAX_IT_GS'
            opt.pf.gs.max_it = ov.(ff);
        case 'ENFORCE_Q_LIMS'
            opt.pf.enforce_q_lims = ov.(ff);
        case 'PF_DC'
            switch ov.(ff)
                case 0
                    opt.model = 'AC';
                case 1
                    opt.model = 'DC';
                otherwise
                    error(errstr, ov.(ff), ff);
            end
        case 'OPF_ALG'
            switch ov.(ff)
                case 0
                    opt.opf.ac.solver = 'DEFAULT';
                case 500
                    opt.opf.ac.solver = 'MINOPF';
                case 520
                    opt.opf.ac.solver = 'FMINCON';
                case {540, 545}
                    opt.opf.ac.solver = 'PDIPM';
                    if ov.(ff) == 545
                        opt.pdipm.step_control = 1;
                    else
                        opt.pdipm.step_control = 0;
                    end
                case 550
                    opt.opf.ac.solver = 'TRALM';
                case {560, 565}
                    opt.opf.ac.solver = 'MIPS';
                    if ov.(ff) == 565
                        opt.mips.step_control = 1;
                    else
                        opt.mips.step_control = 0;
                    end
                case 580
                    opt.opf.ac.solver = 'IPOPT';
                case 600
                    opt.opf.ac.solver = 'KNITRO';
                otherwise
                    error(errstr, ov.(ff), ff);
            end
        case 'OPF_VIOLATION'
            opt.opf.violation = ov.(ff);
        case 'CONSTR_TOL_X'
            opt.fmincon.tol_x = ov.(ff);
            opt.knitro.tol_x = ov.(ff);
        case 'CONSTR_TOL_F'
            opt.fmincon.tol_f = ov.(ff);
            opt.knitro.tol_f = ov.(ff);
        case 'CONSTR_MAX_IT'
            opt.fmincon.max_it = ov.(ff);
        case 'OPF_FLOW_LIM'
            switch ov.(ff)
                case 0
                    opt.opf.flow_lim = 'S';   %% apparent power (MVA)
                case 1
                    opt.opf.flow_lim = 'P';   %% real power (MW)
                case 2
                    opt.opf.flow_lim = 'I';   %% current magnitude (MVA @ 1 p.u. voltage)
                otherwise
                    error(errstr, ov.(ff), ff);
            end
        case 'OPF_IGNORE_ANG_LIM'
            opt.opf.ignore_angle_lim = ov.(ff);
        case 'OPF_ALG_DC'
            switch ov.(ff)
                case 0
                    opt.opf.dc.solver = 'DEFAULT';
                case 100
                    opt.opf.dc.solver = 'BPMPD';
                case {200, 250}
                    opt.opf.dc.solver = 'MIPS';
                    if ov.(ff) == 250
                        opt.mips.step_control = 1;
                    else
                        opt.mips.step_control = 0;
                    end
                case 300
                    opt.opf.dc.solver = 'OT';     %% QUADPROG, LINPROG
                case 400
                    opt.opf.dc.solver = 'IPOPT';
                case 500
                    opt.opf.dc.solver = 'CPLEX';
                case 600
                    opt.opf.dc.solver = 'MOSEK';
                case 700
                    opt.opf.dc.solver = 'GUROBI';
                otherwise
                    error(errstr, ov.(ff), ff);
            end
        case 'VERBOSE'
            opt.verbose = ov.(ff);
        case 'OUT_ALL'
            opt.out.all = ov.(ff);
        case 'OUT_SYS_SUM'
            opt.out.sys_sum = ov.(ff);
        case 'OUT_AREA_SUM'
            opt.out.area_sum = ov.(ff);
        case 'OUT_BUS'
            opt.out.bus = ov.(ff);
        case 'OUT_BRANCH'
            opt.out.branch = ov.(ff);
        case 'OUT_GEN'
            opt.out.gen = ov.(ff);
        case 'OUT_ALL_LIM'
            opt.out.lim.all = ov.(ff);
        case 'OUT_V_LIM'
            opt.out.lim.v = ov.(ff);
        case 'OUT_LINE_LIM'
            opt.out.lim.line = ov.(ff);
        case 'OUT_PG_LIM'
            opt.out.lim.pg = ov.(ff);
        case 'OUT_QG_LIM'
            opt.out.lim.qg = ov.(ff);
        case 'OUT_FORCE'
            opt.out.force = ov.(ff);
        case 'RETURN_RAW_DER'
            opt.opf.return_raw_der = ov.(ff);
        case 'FMC_ALG'
            opt.fmincon.alg = ov.(ff);
        case 'KNITRO_OPT'
            opt.knitro.opt = ov.(ff);
        case 'IPOPT_OPT'
            opt.ipopt.opt = ov.(ff);
        case 'MNS_FEASTOL'
            opt.minopf.feastol = ov.(ff);
        case 'MNS_ROWTOL'
            opt.minopf.rowtol = ov.(ff);
        case 'MNS_XTOL'
            opt.minopf.xtol = ov.(ff);
        case 'MNS_MAJDAMP'
            opt.minopf.majdamp = ov.(ff);
        case 'MNS_MINDAMP'
            opt.minopf.mindamp = ov.(ff);
        case 'MNS_PENALTY_PARM'
            opt.minopf.penalty = ov.(ff);
        case 'MNS_MAJOR_IT'
            opt.minopf.major_it = ov.(ff);
        case 'MNS_MINOR_IT'
            opt.minopf.minor_it = ov.(ff);
        case 'MNS_MAX_IT'
            opt.minopf.max_it = ov.(ff);
        case 'MNS_VERBOSITY'
            opt.minopf.verbosity = ov.(ff);
        case 'MNS_CORE'
            opt.minopf.core = ov.(ff);
        case 'MNS_SUPBASIC_LIM'
            opt.minopf.supbasic_lim = ov.(ff);
        case 'MNS_MULT_PRICE'
            opt.minopf.mult_price = ov.(ff);
        case 'FORCE_PC_EQ_P0'
            opt.sopf.force_Pc_eq_P0 = ov.(ff);
        case 'PDIPM_FEASTOL'
            opt.mips.feastol = ov.(ff);
            opt.pdipm.feastol = ov.(ff);
        case 'PDIPM_GRADTOL'
            opt.mips.gradtol = ov.(ff);
            opt.pdipm.gradtol = ov.(ff);
        case 'PDIPM_COMPTOL'
            opt.mips.comptol = ov.(ff);
            opt.pdipm.comptol = ov.(ff);
        case 'PDIPM_COSTTOL'
            opt.mips.costtol = ov.(ff);
            opt.pdipm.costtol = ov.(ff);
        case 'PDIPM_MAX_IT'
            opt.mips.max_it = ov.(ff);
            opt.pdipm.max_it = ov.(ff);
        case 'SCPDIPM_RED_IT'
            opt.mips.sc.red_it = ov.(ff);
            opt.pdipm.sc.red_it = ov.(ff);
        case 'TRALM_FEASTOL'
            opt.tralm.feastol = ov.(ff);
        case 'TRALM_PRIMETOL'
            opt.tralm.primaltol = ov.(ff);
        case 'TRALM_DUALTOL'
            opt.tralm.dualtol = ov.(ff);
        case 'TRALM_COSTTOL'
            opt.tralm.costtol = ov.(ff);
        case 'TRALM_MAJOR_IT'
            opt.tralm.major_it = ov.(ff);
        case 'TRALM_MINOR_IT'
            opt.tralm.minor_it = ov.(ff);
        case 'SMOOTHING_RATIO'
            opt.pdipm.sc.smooth_ratio = ov.(ff);
            opt.tralm.smooth_ratio = ov.(ff);
        case 'CPLEX_LPMETHOD'
            opt.cplex.lpmethod = ov.(ff);
        case 'CPLEX_QPMETHOD'
            opt.cplex.qpmethod = ov.(ff);
        case 'CPLEX_OPT'
            opt.cplex.opt = ov.(ff);
        case 'MOSEK_LP_ALG'
            opt.mosek.lp_alg = ov.(ff);
        case 'MOSEK_MAX_IT'
            opt.mosek.max_it = ov.(ff);
        case 'MOSEK_GAP_TOL'
            opt.mosek.gap_tol = ov.(ff);
        case 'MOSEK_MAX_TIME'
            opt.mosek.max_time = ov.(ff);
        case 'MOSEK_NUM_THREADS'
            opt.mosek.num_threads = ov.(ff);
        case 'MOSEK_OPT'
            opt.mosek.opt = ov.(ff);
        case 'GRB_METHOD'
            opt.gurobi.method = ov.(ff);
        case 'GRB_TIMELIMIT'
            opt.gurobi.timelimit = ov.(ff);
        case 'GRB_THREADS'
            opt.gurobi.threads = ov.(ff);
        case 'GRB_OPT'
            opt.gurobi.opt = ov.(ff);
        otherwise
            error('mpoption: ''%s'' is not a valid old-style option name', ff);
    end
end
% ov


%%-------------------------------------------------------------------
function opt_s = mpoption_v2s(opt_v)
if DEBUG, fprintf('mpoption_v2s()\n'); end
opt_s = mpoption_default();
errstr = 'mpoption: %g is not a valid value for the old-style ''%s'' option';
switch opt_v(1)                                 %% PF_ALG
    case 1
        opt_s.pf.alg = 'NR';        %% Newton's method
    case 2
        opt_s.pf.alg = 'FDXB';      %% fast-decoupled (XB version)
    case 3
        opt_s.pf.alg = 'FDBX';      %% fast-decoupled (BX version)
    case 4
        opt_s.pf.alg = 'GS';        %% Gauss-Seidel
    otherwise
        error(errstr, opt_v(1), 'PF_ALG');
end
opt_s.pf.tol                = opt_v(2);         %% PF_TOL
opt_s.pf.nr.max_it          = opt_v(3);         %% PF_MAX_IT
opt_s.pf.fd.max_it          = opt_v(4);         %% PF_MAX_IT_FD
opt_s.pf.gs.max_it          = opt_v(5);         %% PF_MAX_IT_GS
opt_s.pf.enforce_q_lims     = opt_v(6);         %% ENFORCE_Q_LIMS
switch opt_v(10)                                %% PF_DC
    case 0
        opt_s.model = 'AC';
    case 1
        opt_s.model = 'DC';
    otherwise
        error(errstr, opt_v(10), 'PF_DC');
end
switch opt_v(11)                                %% OPF_ALG
    case 0
        opt_s.opf.ac.solver = 'DEFAULT';
    case 500
        opt_s.opf.ac.solver = 'MINOPF';
    case 520
        opt_s.opf.ac.solver = 'FMINCON';
    case {540, 545}
        opt_s.opf.ac.solver = 'PDIPM';
    case 550
        opt_s.opf.ac.solver = 'TRALM';
    case {560, 565}
        opt_s.opf.ac.solver = 'MIPS';
    case 580
        opt_s.opf.ac.solver = 'IPOPT';
    case 600
        opt_s.opf.ac.solver = 'KNITRO';
    otherwise
        error(errstr, opt_v(11), 'OPF_ALG');
end
opt_s.opf.violation         = opt_v(16);        %% OPF_VIOLATION

opt_s.fmincon.tol_x         = opt_v(17);        %% CONSTR_TOL_X
opt_s.fmincon.tol_f         = opt_v(18);        %% CONSTR_TOL_F
opt_s.fmincon.max_it        = opt_v(19);        %% CONSTR_MAX_IT

opt_s.knitro.tol_x          = opt_v(17);        %% CONSTR_TOL_X
opt_s.knitro.tol_f          = opt_v(18);        %% CONSTR_TOL_F

switch opt_v(24)                                %% OPF_FLOW_LIM
    case 0
        opt_s.opf.flow_lim = 'S';   %% apparent power (MVA)
    case 1
        opt_s.opf.flow_lim = 'P';   %% real power (MW)
    case 2
        opt_s.opf.flow_lim = 'I';   %% current magnitude (MVA @ 1 p.u. voltage)
    otherwise
        error(errstr, opt_v(10), 'PF_DC');
end

opt_s.opf.ignore_angle_lim  = opt_v(25);        %% OPF_IGNORE_ANG_LIM

switch opt_v(26)                                %% OPF_ALG_DC
    case 0
        opt_s.opf.dc.solver = 'DEFAULT';
    case 100
        opt_s.opf.dc.solver = 'BPMPD';
    case {200, 250}
        opt_s.opf.dc.solver = 'MIPS';
    case 300
        opt_s.opf.dc.solver = 'OT';     %% QUADPROG, LINPROG
    case 400
        opt_s.opf.dc.solver = 'IPOPT';
    case 500
        opt_s.opf.dc.solver = 'CPLEX';
    case 600
        opt_s.opf.dc.solver = 'MOSEK';
    case 700
        opt_s.opf.dc.solver = 'GUROBI';
    otherwise
        error(errstr, opt_v(26), 'OPF_ALG_DC');
end

opt_s.verbose               = opt_v(31);        %% VERBOSE
opt_s.out.all               = opt_v(32);        %% OUT_ALL
opt_s.out.sys_sum           = opt_v(33);        %% OUT_SYS_SUM
opt_s.out.area_sum          = opt_v(34);        %% OUT_AREA_SUM
opt_s.out.bus               = opt_v(35);        %% OUT_BUS
opt_s.out.branch            = opt_v(36);        %% OUT_BRANCH
opt_s.out.gen               = opt_v(37);        %% OUT_GEN
opt_s.out.lim.all           = opt_v(38);        %% OUT_ALL_LIM
opt_s.out.lim.v             = opt_v(39);        %% OUT_V_LIM
opt_s.out.lim.line          = opt_v(40);        %% OUT_LINE_LIM
opt_s.out.lim.pg            = opt_v(41);        %% OUT_PG_LIM
opt_s.out.lim.qg            = opt_v(42);        %% OUT_QG_LIM
opt_s.out.force             = opt_v(44);        %% OUT_FORCE

opt_s.opf.return_raw_der    = opt_v(52);        %% RETURN_RAW_DER

opt_s.fmincon.alg           = opt_v(55);        %% FMC_ALG
opt_s.knitro.opt            = opt_v(58);        %% KNITRO_OPT
opt_s.ipopt.opt             = opt_v(60);        %% IPOPT_OPT

opt_s.minopf.feastol        = opt_v(61);        %% MNS_FEASTOL
opt_s.minopf.rowtol         = opt_v(62);        %% MNS_ROWTOL
opt_s.minopf.xtol           = opt_v(63);        %% MNS_XTOL
opt_s.minopf.majdamp        = opt_v(64);        %% MNS_MAJDAMP
opt_s.minopf.mindamp        = opt_v(65);        %% MNS_MINDAMP
opt_s.minopf.penalty        = opt_v(66);        %% MNS_PENALTY_PARM
opt_s.minopf.major_it       = opt_v(67);        %% MNS_MAJOR_IT
opt_s.minopf.minor_it       = opt_v(68);        %% MNS_MINOR_IT
opt_s.minopf.max_it         = opt_v(69);        %% MNS_MAX_IT
opt_s.minopf.verbosity      = opt_v(70);        %% MNS_VERBOSITY
opt_s.minopf.core           = opt_v(71);        %% MNS_CORE
opt_s.minopf.supbasic_lim   = opt_v(72);        %% MNS_SUPBASIC_LIM
opt_s.minopf.mult_price     = opt_v(73);        %% MNS_MULT_PRICE

opt_s.sopf.force_Pc_eq_P0   = opt_v(80);        %% FORCE_PC_EQ_P0, for c3sopf

if (opt_v(11) == 565 && opt_v(10) == 0) || (opt_v(26) == 250 && opt_v(10) == 1)
    opt_s.mips.step_control = 1;
end
opt_s.mips.feastol          = opt_v(81);        %% PDIPM_FEASTOL
opt_s.mips.gradtol          = opt_v(82);        %% PDIPM_GRADTOL
opt_s.mips.comptol          = opt_v(83);        %% PDIPM_COMPTOL
opt_s.mips.costtol          = opt_v(84);        %% PDIPM_COSTTOL
opt_s.mips.max_it           = opt_v(85);        %% PDIPM_MAX_IT
opt_s.mips.sc.red_it        = opt_v(86);        %% SCPDIPM_RED_IT

opt_s.pdipm.feastol         = opt_v(81);        %% PDIPM_FEASTOL
opt_s.pdipm.gradtol         = opt_v(82);        %% PDIPM_GRADTOL
opt_s.pdipm.comptol         = opt_v(83);        %% PDIPM_COMPTOL
opt_s.pdipm.costtol         = opt_v(84);        %% PDIPM_COSTTOL
opt_s.pdipm.max_it          = opt_v(85);        %% PDIPM_MAX_IT
opt_s.pdipm.sc.red_it       = opt_v(86);        %% SCPDIPM_RED_IT
opt_s.pdipm.sc.smooth_ratio = opt_v(93);        %% SMOOTHING_RATIO
if opt_v(11) == 545 && opt_v(10) == 0
    opt_s.pdipm.step_control = 1;
end

opt_s.tralm.feastol         = opt_v(87);        %% TRALM_FEASTOL
opt_s.tralm.primaltol       = opt_v(88);        %% TRALM_PRIMETOL
opt_s.tralm.dualtol         = opt_v(89);        %% TRALM_DUALTOL
opt_s.tralm.costtol         = opt_v(90);        %% TRALM_COSTTOL
opt_s.tralm.major_it        = opt_v(91);        %% TRALM_MAJOR_IT
opt_s.tralm.minor_it        = opt_v(92);        %% TRALM_MINOR_IT
opt_s.tralm.smooth_ratio    = opt_v(93);        %% SMOOTHING_RATIO

opt_s.cplex.lpmethod        = opt_v(95);        %% CPLEX_LPMETHOD
opt_s.cplex.qpmethod        = opt_v(96);        %% CPLEX_QPMETHOD
opt_s.cplex.opt             = opt_v(97);        %% CPLEX_OPT

opt_s.mosek.lp_alg          = opt_v(111);       %% MOSEK_LP_ALG
opt_s.mosek.max_it          = opt_v(112);       %% MOSEK_MAX_IT
opt_s.mosek.gap_tol         = opt_v(113);       %% MOSEK_GAP_TOL
opt_s.mosek.max_time        = opt_v(114);       %% MOSEK_MAX_TIME
opt_s.mosek.num_threads     = opt_v(115);       %% MOSEK_NUM_THREADS
opt_s.mosek.opt             = opt_v(116);       %% MOSEK_OPT

opt_s.gurobi.method         = opt_v(121);       %% GRB_METHOD
opt_s.gurobi.timelimit      = opt_v(122);       %% GRB_TIMELIMIT
opt_s.gurobi.threads        = opt_v(123);       %% GRB_THREADS
opt_s.gurobi.opt            = opt_v(124);       %% GRB_OPT


%%-------------------------------------------------------------------
function opt_v = mpoption_s2v(opt_s)
if DEBUG, fprintf('mpoption_s2v()\n'); end
%% PF_ALG
old = mpoption_old;
switch upper(opt_s.pf.alg)
    case 'NR'
        PF_ALG = 1;
    case 'FDXB'
        PF_ALG = 2;
    case 'FDBX'
        PF_ALG = 3;
    case 'GS'
        PF_ALG = 4;
end

%% PF_DC
if strcmp(upper(opt_s.model), 'DC')
    PF_DC = 1;
else
    PF_DC = 0;
end

%% OPF_ALG
switch upper(opt_s.opf.ac.solver)
    case 'DEFAULT'
        OPF_ALG = 0;
    case 'MINOPF'
        OPF_ALG = 500;
    case 'FMINCON'
        OPF_ALG = 520;
    case 'PDIPM'
        if isfield(opt_s, 'pdipm') && opt_s.pdipm.step_control
            OPF_ALG = 545;
        else
            OPF_ALG = 540;
        end
    case 'TRALM'
        OPF_ALG = 550;
    case 'MIPS'
        if opt_s.mips.step_control
            OPF_ALG = 565;
        else
            OPF_ALG = 560;
        end
    case 'IPOPT'
        OPF_ALG = 580;
    case 'KNITRO'
        OPF_ALG = 600;
end

%% FMINCON, Knitro tol_x, tol_f, max_it
if strcmp(upper(opt_s.opf.ac.solver), 'KNITRO') && isfield(opt_s, 'knitro')
    CONSTR_TOL_X = opt_s.knitro.tol_x;
    CONSTR_TOL_F = opt_s.knitro.tol_f;
elseif isfield(opt_s, 'fmincon')
    CONSTR_TOL_X  = opt_s.fmincon.tol_x;
    CONSTR_TOL_F  = opt_s.fmincon.tol_f;
else
    CONSTR_TOL_X = old(17);
    CONSTR_TOL_F = old(18);
end
if isfield(opt_s, 'fmincon')
    CONSTR_MAX_IT   = opt_s.fmincon.max_it;
    FMC_ALG         = opt_s.fmincon.alg;
else
    CONSTR_MAX_IT   = old(19);
    FMC_ALG         = old(55);
end

%% OPF_FLOW_LIM
switch upper(opt_s.opf.flow_lim)
    case 'S'
        OPF_FLOW_LIM = 0;
    case {'P', '2'}
        OPF_FLOW_LIM = 1;
    case 'I'
        OPF_FLOW_LIM = 2;
end

%% OPF_ALG_DC
switch upper(opt_s.opf.dc.solver)
    case 'DEFAULT'
        OPF_ALG_DC = 0;
    case 'BPMPD'
        OPF_ALG_DC = 100;
    case 'MIPS'
        if opt_s.mips.step_control
            OPF_ALG_DC = 250;
        else
            OPF_ALG_DC = 200;
        end
    case 'OT'
        OPF_ALG_DC = 300;
    case 'IPOPT'
        OPF_ALG_DC = 400;
    case 'CPLEX'
        OPF_ALG_DC = 500;
    case 'MOSEK'
        OPF_ALG_DC = 600;
    case 'GUROBI'
        OPF_ALG_DC = 700;
end

%% KNITRO_OPT
if isfield(opt_s, 'knitro')
    KNITRO_OPT  = opt_s.knitro.opt;
else
    KNITRO_OPT  = old(58);
end

%% IPOPT_OPT
if isfield(opt_s, 'ipopt')
    IPOPT_OPT  = opt_s.ipopt.opt;
else
    IPOPT_OPT  = old(58);
end

%% MINOPF options
if isfield(opt_s, 'minopf')
    MINOPF_OPTS = [
        opt_s.minopf.feastol;   %% 61 - MNS_FEASTOL
        opt_s.minopf.rowtol;    %% 62 - MNS_ROWTOL
        opt_s.minopf.xtol;      %% 63 - MNS_XTOL
        opt_s.minopf.majdamp;   %% 64 - MNS_MAJDAMP
        opt_s.minopf.mindamp;   %% 65 - MNS_MINDAMP
        opt_s.minopf.penalty;   %% 66 - MNS_PENALTY_PARM
        opt_s.minopf.major_it;  %% 67 - MNS_MAJOR_IT
        opt_s.minopf.minor_it;  %% 68 - MNS_MINOR_IT
        opt_s.minopf.max_it;    %% 69 - MNS_MAX_IT
        opt_s.minopf.verbosity; %% 70 - MNS_VERBOSITY
        opt_s.minopf.core;      %% 71 - MNS_CORE
        opt_s.minopf.supbasic_lim;  %% 72 - MNS_SUPBASIC_LIM
        opt_s.minopf.mult_price;%% 73 - MNS_MULT_PRICE
    ];
else
    MINOPF_OPTS = old(61:73);
end

%% FORCE_PC_EQ_P0
if isfield(opt_s, 'sopf') && isfield(opt_s.sopf, 'force_Pc_eq_P0')
    FORCE_PC_EQ_P0 = opt_s.sopf.force_Pc_eq_P0;
else
    FORCE_PC_EQ_P0 = 0;
end

%% PDIPM options
if isfield(opt_s, 'pdipm')
    PDIPM_OPTS = [
        opt_s.pdipm.feastol;    %% 81 - PDIPM_FEASTOL
        opt_s.pdipm.gradtol;    %% 82 - PDIPM_GRADTOL
        opt_s.pdipm.comptol;    %% 83 - PDIPM_COMPTOL
        opt_s.pdipm.costtol;    %% 84 - PDIPM_COSTTOL
        opt_s.pdipm.max_it;     %% 85 - PDIPM_MAX_IT
        opt_s.pdipm.sc.red_it;  %% 86 - SCPDIPM_RED_IT
    ];
else
    PDIPM_OPTS = old(81:86);
end

%% TRALM options
if isfield(opt_s, 'tralm')
    TRALM_OPTS = [
        opt_s.tralm.feastol;    %% 87 - TRALM_FEASTOL
        opt_s.tralm.primaltol;  %% 88 - TRALM_PRIMETOL
        opt_s.tralm.dualtol;    %% 89 - TRALM_DUALTOL
        opt_s.tralm.costtol;    %% 90 - TRALM_COSTTOL
        opt_s.tralm.major_it;   %% 91 - TRALM_MAJOR_IT
        opt_s.tralm.minor_it;   %% 92 - TRALM_MINOR_IT
    ];
else
    TRALM_OPTS = old(87:92);
end

%% SMOOTHING_RATIO
if strcmp(upper(opt_s.opf.ac.solver), 'TRALM') && isfield(opt_s, 'tralm')
    SMOOTHING_RATIO = opt_s.tralm.smooth_ratio;
elseif isfield(opt_s, 'pdipm')
    SMOOTHING_RATIO = opt_s.pdipm.sc.smooth_ratio;
else
    SMOOTHING_RATIO = old(93);
end

%% CPLEX options
if isfield(opt_s, 'cplex')
    CPLEX_OPTS = [
        opt_s.cplex.lpmethod;   %% 95 - CPLEX_LPMETHOD
        opt_s.cplex.qpmethod;   %% 96 - CPLEX_QPMETHOD
        opt_s.cplex.opt;        %% 97 - CPLEX_OPT
    ];
else
    CPLEX_OPTS = old(95:97);
end

%% MOSEK options
if isfield(opt_s, 'mosek')
    MOSEK_OPTS = [
        opt_s.mosek.lp_alg;     %% 111 - MOSEK_LP_ALG
        opt_s.mosek.max_it;     %% 112 - MOSEK_MAX_IT
        opt_s.mosek.gap_tol;    %% 113 - MOSEK_GAP_TOL
        opt_s.mosek.max_time;   %% 114 - MOSEK_MAX_TIME
        opt_s.mosek.num_threads;%% 115 - MOSEK_NUM_THREADS
        opt_s.mosek.opt;        %% 116 - MOSEK_OPT
    ];
else
    MOSEK_OPTS = old(111:116);
end

%% Gurobi options
if isfield(opt_s, 'gurobi')
    GUROBI_OPTS = [
        opt_s.gurobi.method;    %% 121 - GRB_METHOD
        opt_s.gurobi.timelimit; %% 122 - GRB_TIMELIMIT
        opt_s.gurobi.threads;   %% 123 - GRB_THREADS
        opt_s.gurobi.opt;       %% 124 - GRB_OPT
    ];
else
    GUROBI_OPTS = old(121:124);
end

opt_v = [
        %% power flow options
        PF_ALG;                 %% 1  - PF_ALG
        opt_s.pf.tol;           %% 2  - PF_TOL
        opt_s.pf.nr.max_it;     %% 3  - PF_MAX_IT
        opt_s.pf.fd.max_it;     %% 4  - PF_MAX_IT_FD
        opt_s.pf.gs.max_it;     %% 5  - PF_MAX_IT_GS
        opt_s.pf.enforce_q_lims;%% 6  - ENFORCE_Q_LIMS
        0;                      %% 7  - RESERVED7
        0;                      %% 8  - RESERVED8
        0;                      %% 9  - RESERVED9
        PF_DC;                  %% 10 - PF_DC
        
        %% OPF options
        OPF_ALG;                %% 11 - OPF_ALG
        0;                      %% 12 - RESERVED12 (was OPF_ALG_POLY = 100)
        0;                      %% 13 - RESERVED13 (was OPF_ALG_PWL = 200)
        0;                      %% 14 - RESERVED14 (was OPF_POLY2PWL_PTS = 10)
        0;                      %% 15 - OPF_NEQ (removed)
        opt_s.opf.violation;    %% 16 - OPF_VIOLATION
        CONSTR_TOL_X;           %% 17 - CONSTR_TOL_X
        CONSTR_TOL_F;           %% 18 - CONSTR_TOL_F
        CONSTR_MAX_IT;          %% 19 - CONSTR_MAX_IT
        old(20);                %% 20 - LPC_TOL_GRAD (removed)
        old(21);                %% 21 - LPC_TOL_X (removed)
        old(22);                %% 22 - LPC_MAX_IT (removed)
        old(23);                %% 23 - LPC_MAX_RESTART (removed)
        OPF_FLOW_LIM;           %% 24 - OPF_FLOW_LIM
        opt_s.opf.ignore_angle_lim; %% 25 - OPF_IGNORE_ANG_LIM
        OPF_ALG_DC;             %% 26 - OPF_ALG_DC
        0;                      %% 27 - RESERVED27
        0;                      %% 28 - RESERVED28
        0;                      %% 29 - RESERVED29
        0;                      %% 30 - RESERVED30
        
        %% output options
        opt_s.verbose;          %% 31 - VERBOSE
        opt_s.out.all;          %% 32 - OUT_ALL
        opt_s.out.sys_sum;      %% 33 - OUT_SYS_SUM
        opt_s.out.area_sum;     %% 34 - OUT_AREA_SUM
        opt_s.out.bus;          %% 35 - OUT_BUS
        opt_s.out.branch;       %% 36 - OUT_BRANCH
        opt_s.out.gen;          %% 37 - OUT_GEN
        opt_s.out.lim.all;      %% 38 - OUT_ALL_LIM
        opt_s.out.lim.v;        %% 39 - OUT_V_LIM
        opt_s.out.lim.line;     %% 40 - OUT_LINE_LIM
        opt_s.out.lim.pg;       %% 41 - OUT_PG_LIM
        opt_s.out.lim.qg;       %% 42 - OUT_QG_LIM
        0;                      %% 43 - RESERVED43 (was OUT_RAW)
        opt_s.out.force;        %% 44 - OUT_FORCE
        0;                      %% 45 - RESERVED45
        0;                      %% 46 - RESERVED46
        0;                      %% 47 - RESERVED47
        0;                      %% 48 - RESERVED48
        0;                      %% 49 - RESERVED49
        0;                      %% 50 - RESERVED50
        
        %% other options
        old(51);                %% 51 - SPARSE_QP (removed)
        opt_s.opf.return_raw_der;   %% 52 - RETURN_RAW_DER
        0;                      %% 53 - RESERVED53
        0;                      %% 54 - RESERVED54
        FMC_ALG;                %% 55 - FMC_ALG
        0;                      %% 56 - RESERVED56
        0;                      %% 57 - RESERVED57
        KNITRO_OPT;             %% 58 - KNITRO_OPT
        0;                      %% 59 - RESERVED59
        IPOPT_OPT;              %% 60 - IPOPT_OPT
        
        %% MINOPF options
        MINOPF_OPTS;            %% 61-73 - MNS_FEASTOL-MNS_MULT_PRICE
        0;                      %% 74 - RESERVED74
        0;                      %% 75 - RESERVED75
        0;                      %% 76 - RESERVED76
        0;                      %% 77 - RESERVED77
        0;                      %% 78 - RESERVED78
        0;                      %% 79 - RESERVED79
        FORCE_PC_EQ_P0;         %% 80 - FORCE_PC_EQ_P0, for c3sopf
        
        %% MIPS, PDIPM, SC-PDIPM, and TRALM options
        PDIPM_OPTS;             %% 81-86 - PDIPM_FEASTOL-SCPDIPM_RED_IT
        TRALM_OPTS;             %% 87-92 - TRALM_FEASTOL-TRALM_MINOR_IT
        SMOOTHING_RATIO;        %% 93 - SMOOTHING_RATIO
        0;                      %% 94 - RESERVED94
        
        %% CPLEX options
        CPLEX_OPTS;             %% 95-97 - CPLEX_LPMETHOD-CPLEX_OPT
        0;                      %% 98 - RESERVED98
        0;                      %% 99 - RESERVED99
        0;                      %% 100 - RESERVED100
        0;                      %% 101 - RESERVED101
        0;                      %% 102 - RESERVED102
        0;                      %% 103 - RESERVED103
        0;                      %% 104 - RESERVED104
        0;                      %% 105 - RESERVED105
        0;                      %% 106 - RESERVED106
        0;                      %% 107 - RESERVED107
        0;                      %% 108 - RESERVED108
        0;                      %% 109 - RESERVED109
        0;                      %% 110 - RESERVED110

        %% MOSEK options
        MOSEK_OPTS;             %% 111-116 - MOSEK_LP_ALG-MOSEK_OPT
        0;                      %% 117 - RESERVED117
        0;                      %% 118 - RESERVED118
        0;                      %% 119 - RESERVED119
        0;                      %% 120 - RESERVED120

        %% Gurobi options
        GUROBI_OPTS;            %% 121-124 - GRB_METHOD-GRB_OPT
    ];


%%-------------------------------------------------------------------
function optt = mpoption_default()
if DEBUG, fprintf('mpoption_default()\n'); end
persistent opt;             %% cache this for speed
if ~isstruct(opt)
    opt = struct(...
        'v',                    mpoption_version, ...   %% version
        'model',                'AC', ...
        'pf',                   struct(...
            'alg',                  'NR', ...
            'tol',                  1e-8, ...
            'nr',                   struct(...
                'max_it',               10, ...
                'lin_solver',            '' ), ...
            'fd',                   struct(...
                'max_it',               30  ), ...
            'gs',                   struct(...
                'max_it',               1000    ), ...
            'radial',               struct(...
                'max_it',               20   , ...
                'vcorr',                 0  ), ...
            'enforce_q_lims',       0   ), ...
        'cpf',                  struct(...
            'parameterization',     3, ...
            'stop_at',              'NOSE', ...     %% 'NOSE', <lam val>, 'FULL'
            'enforce_p_lims',       0, ...
            'enforce_q_lims',       0, ...
            'enforce_v_lims',       0, ...
            'enforce_flow_lims',    0, ...
            'step',                 0.05, ...
            'step_min',             1e-4, ...
            'step_max',             0.2, ...
            'adapt_step',           0, ...
            'adapt_step_damping',   0.7, ...
            'adapt_step_tol',       1e-3, ...
            'target_lam_tol',       1e-5, ...
            'nose_tol',             1e-5, ...
            'p_lims_tol',           0.01, ...
            'q_lims_tol',           0.01, ...
            'v_lims_tol',           1e-4, ...
            'flow_lims_tol',        0.01, ...
            'plot',                 struct(...
                'level',                0, ...
                'bus',                  []  ), ...
            'user_callback',        ''  ), ...
        'opf',                  struct(...
            'ac',                   struct(...
                'solver',               'DEFAULT'   ), ...
            'dc',                   struct(...
                'solver',               'DEFAULT'   ), ...
            'current_balance',      0, ...
            'v_cartesian',          0, ...
            'violation',            5e-6, ...
            'use_vg',               0, ...
            'flow_lim',             'S', ...
            'ignore_angle_lim',     0, ...
            'softlims',             struct(...
                'default',               1  ), ...
            'init_from_mpc',        -1, ...
            'start',                0, ...
            'return_raw_der',       0   ), ...
        'verbose',              1, ...
        'out',                  struct(...
            'all',                  -1, ...
            'sys_sum',              1, ...
            'area_sum',             0, ...
            'bus',                  1, ...
            'branch',               1, ...
            'gen',                  0, ...
            'lim',                  struct(...
                'all',                  -1, ...
                'v',                    1, ...
                'line',                 1, ...
                'pg',                   1, ...
                'qg',                   1   ), ...
            'force',                0, ...
            'suppress_detail',      -1  ), ...
        'mips',                 struct(...  %% see mpoption_info_mips() for optional fields
            'step_control',         0, ...
            'linsolver',            '', ...
            'feastol',              0, ...
            'gradtol',              1e-6, ...
            'comptol',              1e-6, ...
            'costtol',              1e-6, ...
            'max_it',               150, ...
            'sc',                   struct(...
                'red_it',               20  )), ...
        'exp',                  struct(... %% experimental options
            'sys_wide_zip_loads',   struct(...
                'pw',                   [], ...
                'qw',                   []  )) ...
    );
    opt_pkgs = mpoption_optional_pkgs();
    for k = 1:length(opt_pkgs)
        fname = ['mpoption_info_' opt_pkgs{k}];
        if exist(fname, 'file') == 2
            opt = nested_struct_copy(opt, feval(fname, 'D'));
        end
    end
end
optt = opt;

%%-------------------------------------------------------------------
function optt = mpoption_optional_fields()
if DEBUG, fprintf('mpoption_optional_fields()\n'); end
persistent opt;         %% cache this for speed
if ~isstruct(opt)
    opt_pkgs = mpoption_optional_pkgs();
    opt = struct;
    for k = 1:length(opt_pkgs)
        fname = ['mpoption_info_' opt_pkgs{k}];
        if exist(fname, 'file') == 2
            opt = nested_struct_copy(opt, feval(fname, 'V'));
        end
    end
end
optt = opt;

%% globals
%%-------------------------------------------------------------------
function v = mpoption_version
v = 19;     %% version number of MATPOWER options struct
            %% (must be incremented every time structure is updated)
            %% v1   - first version based on struct (MATPOWER 5.0b1)
            %% v2   - added 'linprog' and 'quadprog' fields
            %% v3   - (forgot to increment v) added 'out.suppress_detail'
            %%        field
            %% v4   - (forgot to increment v) MIPS 1.1, added optional
            %%        fields to 'mips' options: xi, sigma, z0, alpha_min,
            %%        rho_min, rho_max, mu_threshold and max_stepsize
            %% v5   - (forgot to increment v) added 'opf.init_from_mpc'
            %%        field (MATPOWER 5.0)
            %% v6   - added 'clp' field
            %% v7   - added 'intlinprog' field
            %% v8   - MIPS 1.2, added 'linsolver' field to
            %%        'mips' options
            %% v9   - added 'exp' for experimental fields, specifically
            %%        'sys_wide_zip_loads.pw', 'sys_wide_zip_loads.qw'
            %% v10  - added 'most' field
            %% v11  - added 'cpf' options 'adapt_step_damping',
            %%        'enforce_p_lims', 'enforce_q_lims', 'target_lam_tol'
            %%        'nose_tol', 'p_lims_tol' and 'q_lims_tol',
            %%        removed option 'cpf.user_callback_args'
            %% v12  - added option 'opf.use_vg'
            %% v13  - added 'pf.radial.max_it', 'pf.radial.vcorr'
            %% v14  - added 'pf.nr.lin_solver'
            %% v15  - added 'cpf.enforce_v_lims', 'cpf.enforce_flow_lims',
            %%        'cpf.v_lims_tol', and 'cpf.flow_lims_tol'
            %% v16  - added 'opf.start' (deprecated 'opf.init_from_mpc')
            %% v17  - added 'knitro.maxit'
            %% v18  - added 'opf.current_balance' and 'opf.v_cartesian'
            %% v19  - added 'opf.softlims.default'

%%-------------------------------------------------------------------
function db_level = DEBUG
db_level = 0;

%%-------------------------------------------------------------------
function pkgs = mpoption_optional_pkgs()
pkgs = {...
    'clp', 'cplex', 'fmincon', 'gurobi', 'glpk', 'intlinprog', 'ipopt', ...
    'knitro', 'linprog', 'minopf', 'most', 'mosek', 'quadprog', 'sdp_pf', ...
    'sopf', 'tspopf', 'yalmip' ...
};
