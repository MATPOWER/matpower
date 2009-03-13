function [options, names] = mpoption(varargin)
%MPOPTION  Used to set and retrieve a MATPOWER options vector.
%
%   opt = mpoption
%       returns the default options vector
%
%   opt = mpoption(name1, value1, name2, value2, ...)
%       returns the default options vector with new values for up to 7
%       options, name# is the name of an option, and value# is the new
%       value. Example: options = mpoption('PF_ALG', 2, 'PF_TOL', 1e-4)
%
%   opt = mpoption(opt, name1, value1, name2, value2, ...)
%       same as above except it uses the options vector opt as a base
%       instead of the default options vector.
%
%   The currently defined options are as follows:
%
%      idx - NAME, default          description [options]
%      ---   -------------          -------------------------------------
%   power flow options
%       1  - PF_ALG, 1              power flow algorithm
%           [   1 - Newton's method                                     ]
%           [   2 - Fast-Decoupled (XB version)                         ]
%           [   3 - Fast-Decoupled (BX version)                         ]
%           [   4 - Gauss Seidel                                        ]
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
%           [    0 - do NOT enforce limits                              ]
%           [    1 - enforce limits, simultaneous bus type conversion   ]
%           [    2 - enforce limits, one-at-a-time bus type conversion  ]
%       10 - PF_DC, 0               use DC power flow formulation, for
%                                   power flow and OPF
%           [    0 - use AC formulation & corresponding algorithm opts  ]
%           [    1 - use DC formulation, ignore AC algorithm options    ]
%   OPF options
%       11 - OPF_ALG, 0             algorithm to use for OPF
%           [    0 - choose best default solver available in the        ]
%           [        following order, 500, 540, 560                     ]
%           [  300 - generalized formulation, constr                    ]
%           [  320 - generalized formulation, dense LP                  ]
%           [  340 - generalized formulation, sparse LP (relaxed)       ]
%           [  360 - generalized formulation, sparse LP (full)          ]
%           [  500 - generalized formulation, MINOS                     ]
%           [  520 - generalized formulation, fmincon                   ]
%           [  540 - generalized formulation, PDIPM                     ]
%           [        primal/dual interior point method                  ]
%           [  545 - generalized formulation (except CCV), SCPDIPM      ]
%           [        step-controlled primal/dual interior point method  ]
%           [  550 - generalized formulation, TRALM                     ]
%           [        trust region based augmented Langrangian method    ]
%           [  560 - generalized formulation, PDIPM (pure Matlab)       ]
%           [        primal/dual interior point method                  ]
%           [ See the User's Manual for details on the formulations.    ]
%       16 - OPF_VIOLATION, 5e-6    constraint violation tolerance
%       17 - CONSTR_TOL_X, 1e-4     termination tol on x for copf & fmincopf
%       18 - CONSTR_TOL_F, 1e-4     termination tol on F for copf & fmincopf
%       19 - CONSTR_MAX_IT, 0       max number of iterations for copf & fmincopf
%                                   [       0 => 2*nb + 150             ]
%       20 - LPC_TOL_GRAD, 3e-3     termination tolerance on gradient for lpopf
%       21 - LPC_TOL_X, 1e-4        termination tolerance on x (min step size)
%                                   for lpopf
%       22 - LPC_MAX_IT, 400        maximum number of iterations for lpopf
%       23 - LPC_MAX_RESTART, 5     maximum number of restarts for lpopf
%       24 - OPF_FLOW_LIM, 0        qty to limit for branch flow constraints
%           [   0 - apparent power flow (limit in MVA)                  ]
%           [   1 - active power flow (limit in MW)                     ]
%           [   2 - current magnitude (limit in MVA at 1 p.u. voltage   ]
%       25 - OPF_IGNORE_ANG_LIM, 0  ignore angle difference limits for branches
%                                   even if specified       [   0 or 1  ]
%   output options
%       31 - VERBOSE, 1             amount of progress info printed
%           [   0 - print no progress info                              ]
%           [   1 - print a little progress info                        ]
%           [   2 - print a lot of progress info                        ]
%           [   3 - print all progress info                             ]
%       32 - OUT_ALL, -1            controls printing of results
%           [  -1 - individual flags control what prints                ]
%           [   0 - don't print anything                                ]
%           [       (overrides individual flags, except OUT_RAW)        ]
%           [   1 - print everything                                    ]
%           [       (overrides individual flags, except OUT_RAW)        ]
%       33 - OUT_SYS_SUM, 1         print system summary    [   0 or 1  ]
%       34 - OUT_AREA_SUM, 0        print area summaries    [   0 or 1  ]
%       35 - OUT_BUS, 1             print bus detail        [   0 or 1  ]
%       36 - OUT_BRANCH, 1          print branch detail     [   0 or 1  ]
%       37 - OUT_GEN, 0             print generator detail  [   0 or 1  ]
%                                   (OUT_BUS also includes gen info)
%       38 - OUT_ALL_LIM, -1        control constraint info output
%           [  -1 - individual flags control what constraint info prints]
%           [   0 - no constraint info (overrides individual flags)     ]
%           [   1 - binding constraint info (overrides individual flags)]
%           [   2 - all constraint info (overrides individual flags)    ]
%       39 - OUT_V_LIM, 1           control output of voltage limit info
%           [   0 - don't print                                         ]
%           [   1 - print binding constraints only                      ]
%           [   2 - print all constraints                               ]
%           [   (same options for OUT_LINE_LIM, OUT_PG_LIM, OUT_QG_LIM) ]
%       40 - OUT_LINE_LIM, 1        control output of line limit info
%       41 - OUT_PG_LIM, 1          control output of gen P limit info
%       42 - OUT_QG_LIM, 1          control output of gen Q limit info
%       43 - OUT_RAW, 0             print raw data for Perl database
%                                   interface code          [   0 or 1  ]
%   other options
%       51 - SPARSE_QP, 1           pass sparse matrices to QP and LP
%                                   solvers if possible     [   0 or 1  ]
%   fmincon options
%       55 - FMC_ALG, 1             algorithm used by fmincon for OPF
%                                   for Optimization Toolbox 4 and later
%            [  1 - active-set                                          ]
%            [  2 - interior-point, w/default 'bfgs' Hessian approx     ]
%            [  3 - interior-point, w/ 'lbfgs' Hessian approx           ]
%            [  4 - interior-point, w/exact user-supplied Hessian       ]
%            [  5 - interior-point, w/Hessian via finite differences    ]
%
%   MINOPF options
%       61 - MNS_FEASTOL, 0 (1E-3)  primal feasibility tolerance,
%                                   set to value of OPF_VIOLATION by default
%       62 - MNS_ROWTOL, 0  (1E-3)  row tolerance
%                                   set to value of OPF_VIOLATION by default
%       63 - MNS_XTOL, 0     (1E-3) x tolerance
%                                   set to value of CONSTR_TOL_X by default
%       64 - MNS_MAJDAMP, 0 (0.5)   major damping parameter
%       65 - MNS_MINDAMP, 0 (2.0)   minor damping parameter
%       66 - MNS_PENALTY_PARM, 0 (1.0)  penalty parameter
%       67 - MNS_MAJOR_IT, 0 (200)  major iterations
%       68 - MNS_MINOR_IT, 0 (2500) minor iterations
%       69 - MNS_MAX_IT, 0 (2500)   iterations limit
%       70 - MNS_VERBOSITY, -1
%           [  -1 - controlled by VERBOSE flag                          ]
%           [   0 - print nothing                                       ]
%           [   1 - print only termination status message               ]
%           [   2 - print termination status and screen progress        ]
%           [   3 - print screen progress, report file (usually fort.9) ]
%       71 - MNS_CORE, 1200 * nb + 2 * (nb + ng)^2 
%       72 - MNS_SUPBASIC_LIM, 0 (2*nb + 2*ng) superbasics limit
%       73 - MNS_MULT_PRICE, 0 (30) multiple price
%
%   PDIPM, SC-PDIPM, and TRALM options
%       81 - PDIPM_FEASTOL, 0       feasibility (equality) tolerance for 
%                                   PDIPM and SC-PDIPM
%                                   set to value of OPF_VIOLATION by default
%       82 - PDIPM_GRADTOL, 1e-6    gradient tolerance for PDIPM 
%                                   and SC-PDIPM
%       83 - PDIPM_COMPTOL, 1e-6    complementary condition (inequality) 
%                                   tolerance for PDIPM and SC-PDIPM
%       84 - PDIPM_COSTTOL, 1e-6    optimality tolerance for PDIPM and 
%                                   SC-PDIPM
%       85 - PDIPM_MAX_IT,  150     maximum number of iterations for 
%                                   PDIPM and SC-PDIPM
%       86 - SCPDIPM_RED_IT, 20     maximum number of SC-PDIPM reductions 
%                                   per iteration
%       87 - TRALM_FEASTOL, 0       feasibility tolerance for TRALM
%                                   set to value of OPF_VIOLATION by default
%       88 - TRALM_PRIMETOL, 5e-4   prime variable tolerance for TRALM
%       89 - TRALM_DUALTOL, 5e-4    dual variable tolerance for TRALM
%       90 - TRALM_COSTTOL, 1e-5    optimality tolerance for TRALM
%       91 - TRALM_MAJOR_IT, 40     maximum number of major iterations
%       92 - TRALM_MINOR_IT, 100    maximum number of minor iterations
%       93 - SMOOTHING_RATIO, 0.04  piecewise linear curve smoothing ratio
%                                   used in SC-PDIPM and TRALM

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
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
        100;    %% 12 - OPF_ALG_POLY, deprecated
        200;    %% 13 - OPF_ALG_PWL, deprecated
        10;     %% 14 - OPF_POLY2PWL_PTS, deprecated
        0;      %% 15 - OPF_NEQ, not a user option (number of eq constraints for
                %%          copf, lpopf and dcopf algorithms, set by program)
        5e-6;   %% 16 - OPF_VIOLATION
        1e-4;   %% 17 - CONSTR_TOL_X
        1e-4;   %% 18 - CONSTR_TOL_F
        0;      %% 19 - CONSTR_MAX_IT
        3e-3;   %% 20 - LPC_TOL_GRAD
        1e-4;   %% 21 - LPC_TOL_X
        400;    %% 22 - LPC_MAX_IT
        5;      %% 23 - LPC_MAX_RESTART
        0;      %% 24 - OPF_FLOW_LIM
        0;      %% 25 - OPF_IGNORE_ANG_LIM
        0;      %% 26 - RESERVED26
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
        0;      %% 43 - OUT_RAW
        0;      %% 44 - RESERVED44
        0;      %% 45 - RESERVED45
        0;      %% 46 - RESERVED46
        0;      %% 47 - RESERVED47
        0;      %% 48 - RESERVED48
        0;      %% 49 - RESERVED49
        0;      %% 50 - RESERVED50
        
        %% other options
        1;      %% 51 - SPARSE_QP
        0;      %% 52 - RESERVED52
        0;      %% 53 - RESERVED53
        0;      %% 54 - RESERVED54
        1;      %% 55 - FMC_ALG
        0;      %% 56 - RESERVED56
        0;      %% 57 - RESERVED57
        0;      %% 58 - RESERVED58
        0;      %% 59 - RESERVED59
        0;      %% 60 - RESERVED60
        
        %% other options
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
        
        %% PDIPM options
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
    ];
end

%%-----  set up option names  -----
%% power flow options
names = str2mat(    'PF_ALG', ...               %% 1
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
names = str2mat(    names, ...
                    'OPF_ALG', ...              %% 11
                    'OPF_ALG_POLY', ...         %% 12
                    'OPF_ALG_PWL', ...          %% 13
                    'OPF_POLY2PWL_PTS', ...     %% 14
                    'OPF_NEQ', ...              %% 15
                    'OPF_VIOLATION', ...        %% 16
                    'CONSTR_TOL_X', ...         %% 17
                    'CONSTR_TOL_F', ...         %% 18
                    'CONSTR_MAX_IT', ...        %% 19
                    'LPC_TOL_GRAD'  );          %% 20
names = str2mat(    names, ...
                    'LPC_TOL_X', ...            %% 21
                    'LPC_MAX_IT', ...           %% 22
                    'LPC_MAX_RESTART', ...      %% 23
                    'OPF_FLOW_LIM', ...         %% 24
                    'OPF_IGNORE_ANG_LIM', ...   %% 25
                    'RESERVED26', ...           %% 26
                    'RESERVED27', ...           %% 27
                    'RESERVED28', ...           %% 28
                    'RESERVED29', ...           %% 29
                    'RESERVED30'    );          %% 30

%% output options
names = str2mat(    names, ...
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
names = str2mat(    names, ...
                    'OUT_PG_LIM', ...           %% 41
                    'OUT_QG_LIM', ...           %% 42
                    'OUT_RAW', ...              %% 43
                    'RESERVED44', ...           %% 44
                    'RESERVED45', ...           %% 45
                    'RESERVED46', ...           %% 46
                    'RESERVED47', ...           %% 47
                    'RESERVED48', ...           %% 48
                    'RESERVED49', ...           %% 49
                    'RESERVED50'    );          %% 50
%% other options
names = str2mat(    names, ...
                    'SPARSE_QP', ...            %% 51
                    'RESERVED52', ...           %% 52
                    'RESERVED53', ...           %% 53
                    'RESERVED54', ...           %% 54
                    'FMC_ALG', ...              %% 55
                    'RESERVED56', ...           %% 56
                    'RESERVED57', ...           %% 57
                    'RESERVED58', ...           %% 58
                    'RESERVED59', ...           %% 59
                    'RESERVED60'    );          %% 60
%% MINOS options
names = str2mat(    names, ...
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
names = str2mat(    names, ...
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

%% PDIPM, SC-PDIPM, and TRALM options                
names = str2mat(    names, ...
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
                    'SMOOTHING_RATIO'    );     %% 93
                
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
        error(sprintf('"%s" is not a valid named option', pname));
    end
    % fprintf('''%s'' (%d) = %d\n', pname, pidx, pval);

    %% update option
    options(pidx) = pval;

    i = i + 2;                              %% go to next parameter
end

return;
