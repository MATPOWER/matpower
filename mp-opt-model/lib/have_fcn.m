function varargout = have_fcn(varargin)
%HAVE_FCN  Test for optional functionality / version info.
%
%   -----  DEPRECATED - Please use HAVE_FEATURE instead    -----
%%
%   TORF = HAVE_FCN(TAG)
%   TORF = HAVE_FCN(TAG, TOGGLE)
%   VER_STR = HAVE_FCN(TAG, 'vstr')
%   VER_NUM = HAVE_FCN(TAG, 'vnum')
%   DATE    = HAVE_FCN(TAG, 'date')
%   INFO    = HAVE_FCN(TAG, 'all')
%   HAVE_FCN(TAG, 'clear_cache')
%   HAVE_FCN('all', 'clear_cache')
%
%   Returns availability, version and release information for optional
%   functionality. All information is cached, and the cached values
%   returned on subsequent calls. If the functionality exists, an attempt is
%   made to determine the release date and version number. The second
%   argument defines which value is returned, as follows:
%       <none>      1 = optional functionality is available, 0 = not available
%       'vstr'      version number as a string (e.g. '3.11.4')
%       'vnum'      version number as numeric value (e.g. 3.011004)
%       'date'      release date as a string (e.g. '20-Jan-2015')
%       'all'       struct with fields named 'av' (for 'availability'), 'vstr',
%                   'vnum' and 'date', and values corresponding to the above,
%                   respectively.
%
%   For functionality that is not available, all calls with a string-valued
%   second argument will return an empty value.
%
%   Alternatively, the optional functionality specified by TAG can be toggled
%   OFF or ON by calling HAVE_FCN with a numeric second argument TOGGLE with
%   one of the following values:
%       0 - turn OFF the optional functionality
%       1 - turn ON the optional functionality (if available)
%      -1 - toggle the ON/OFF state of the optional functionality
%
%   Finally, passing 'clear_cache' as the second argument will cause the
%   cached information to be cleared for the specified TAG or, if the first
%   argument is 'all', for all optional functionality. When calling with
%   'clear_cache' no return value is defined.
%
%   Possible values for input TAG and their meanings:
%       bpmpd       - BP, BPMPD interior point LP/QP solver
%       clp         - CLP, LP/QP solver(https://github.com/coin-or/Clp)
%         opti_clp  - version of CLP distributed with OPTI Toolbox
%                       (https://www.inverseproblem.co.nz/OPTI/)
%       cplex       - CPLEX, IBM ILOG CPLEX Optimizer
%       fmincon     - FMINCON, solver from Optimization Toolbox
%         fmincon_ipm - FMINCON with Interior Point solver, from Opt Tbx 4.x+
%       fsolve      - FSOLVE, nonlinear equation solver from Opt Toolbox
%       glpk        - GLPK, GNU Linear Programming Kit
%       gurobi      - GUROBI, Gurobi solver (https://www.gurobi.com/)
%       intlinprog  - INTLINPROG, MILP solver from Optimization
%                     Toolbox 7.0 (R2014a)+
%       ipopt       - IPOPT, NLP solver
%                       (https://github.com/coin-or/Ipopt)
%       knitro      - Artelys Knitro, NLP solver
%                     (https://www.artelys.com/solvers/knitro/)
%         knitromatlab - Artelys Knitro, version 9.0.0+
%         ktrlink      - KNITRO, version < 9.0.0 (requires Opt Tbx)
%       linprog     - LINPROG, LP solver from Optimization Toolbox
%         linprog_ds - LINPROG with dual-simplex solver
%                       from Optimization Toolbox 7.1 (R2014b) +
%       matlab      - code is running under MATLAB, as opposed to Octave
%       mosek       - MOSEK, LP/QP solver (https://www.mosek.com/)
%       octave      - code is running under GNU Octave, as opposed to MATLAB
%       optim       - Optimization Toolbox
%       optimoptions - OPTIMOPTIONS, option setting funciton for Optim Tbx 6.3+
%       osqp        - OSQP (Operator Splitting QP) solver (https://osqp.org)
%       pardiso     - PARDISO, Parallel Sparse Direct & Iterative Linear Solver
%                       (https://pardiso-project.org)
%       quadprog    - QUADPROG, QP solver from Optimization Toolbox 2.x +
%         quadprog_ls - QUADPROG with large-scale interior point convex solver
%                       from Optimization Toolbox 6.x +
%       sdpt3       - SDPT3 SDP solver (https://github.com/sqlp/sdpt3)
%       sedumi      - SeDuMi SDP solver (http://sedumi.ie.lehigh.edu)
%       yalmip      - YALMIP SDP modeling platform (https://yalmip.github.io)
%
%     Functionality related to MATPOWER
%       e4st        - E4ST (https://e4st.com/)
%       minopf      - MINOPF, MINOPF, MINOS-based OPF solver
%       most        - MOST, MATPOWER Optimal Scheduling Tool
%       pdipmopf    - PDIPMOPF, primal-dual interior point method OPF solver
%       scpdipmopf  - SCPDIPMOPF, step-controlled PDIPM OPF solver
%       sdp_pf      - SDP_PF applications of semi-definite programming
%                     relaxation of power flow equations
%       smartmarket - RUNMARKET and friends, for running an energy auction
%       syngrid     - SynGrid, Synthetic Grid Creation for MATPOWER
%       tralmopf    - TRALMOPF, trust region based augmented Langrangian
%                     OPF solver
%
%   Examples:
%       if have_fcn('minopf')
%           results = runopf(mpc, mpoption('opf.ac.solver', 'MINOPF'));
%       end

%   Private tags for internal use only:
%       catchme         - support for 'catch me' syntax in try/catch constructs
%       evalc           - support for evalc() function
%       ipopt_auxdata   - support for ipopt_auxdata(), required by 3.11 and later
%       lu_vec          - support for lu(..., 'vector') syntax
%       pardiso_legacy  - PARDISO v5, individual MEX files for factor, solve, etc
%       pardiso_object  - PARDISO v6 and later, object interface
%       regexp_split    - support for 'split' argument to regexp()
%       rithmaticker    - used for testing HAVE_FCN
%
%   The following calling syntaxes are also implemented to set and get the
%   entire cache struct and are used during testing only.
%       CACHE = HAVE_FCN('all', 'get_cache')
%       HAVE_FCN(CACHE, 'set_cache')

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

[varargout{1:nargout}] = have_feature(varargin{:});
