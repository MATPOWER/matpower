function [options, names] = mpoption(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15)
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
%       2  - PF_TOL, 1e-8           termination tolerance on per unit
%                                   P & Q mismatch
%       3  - PF_MAX_IT, 10          maximum number of iterations for
%                                   Newton's method
%       4  - PF_MAX_IT_FD, 30       maximum number of iterations for 
%                                   fast decoupled method
%   OPF options
%       11 - OPF_ALG, 0             algorithm to use for OPF
%           [   (see README for more info on formulations/algorithms)   ]
%           [   The algorithm code = F * 100 + S * 20, where ...        ]
%           [       F specifies one of the following OPF formulations   ]
%           [           1 - standard (polynomial cost in obj fcn)       ]
%           [           2 - CCV (constrained cost variables)            ]
%           [       S specifies one of the following solvers            ]
%           [           0 - 'constr' from Optimization Toolbox          ]
%           [           1 - Dense LP-based method                       ]
%           [           2 - Sparse LP-based method w/relaxed constraints]
%           [           3 - Sparse LP-based method w/full constraints   ]
%           [   This yields the following 9 codes:                      ]
%           [         0 - choose appropriate default from OPF_ALG_POLY  ]
%           [             or OPF_ALG_PWL                                ]
%           [       100 - standard formulation, constr                  ]
%           [       120 - standard formulation, dense LP                ]
%           [       140 - standard formulation, sparse LP (relaxed)     ]
%           [       160 - standard formulation, sparse LP (full)        ]
%           [       200 - CCV formulation, constr                       ]
%           [       220 - CCV formulation, dense LP                     ]
%           [       240 - CCV formulation, sparse LP (relaxed)          ]
%           [       260 - CCV formulation, sparse LP (full)             ]
%       12 - OPF_ALG_POLY, 100      default OPF algorithm for use with
%                                   polynomial cost functions
%       13 - OPF_ALG_PWL, 200       default OPF algorithm for use with
%                                   piece-wise linear cost functions
%       14 - OPF_POLY2PWL_PTS, 10   number of evaluation points to use
%                                   when converting from polynomial to
%                                   piece-wise linear costs
%       15 - OPF_NEQ, 0             number of equality constraints
%                                   (0 => 2*nb, set by program, not a
%                                   user option)
%       16 - OPF_VIOLATION, 5e-6    constraint violation tolerance
%       17 - CONSTR_TOL_X, 1e-4     termination tol on x for 'constr'
%       18 - CONSTR_TOL_F, 1e-4     termination tol on F for 'constr'
%       19 - CONSTR_MAX_IT, 0       max number of iterations for 'constr'
%                                   [       0 => 2*nb + 150             ]
%       20 - LPC_TOL_GRAD, 3e-3     termination tolerance on gradient
%                                   for 'LPconstr'
%       21 - LPC_TOL_X, 5e-3        termination tolerance on x (min
%                                   step size) for 'LPconstr'
%       22 - LPC_MAX_IT, 1000       maximum number of iterations for
%                                   'LPconstr'
%       23 - LPC_MAX_RESTART, 5     maximum number of restarts for
%                                   'LPconstr'
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

%   MATPOWER Version 2.0
%   by Ray Zimmerman, PSERC Cornell    12/24/97
%   Copyright (c) 1996, 1997 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%%-----  set up default option values  -----
i = 1;
if rem(nargin, 2)		%% odd number of arguments
	options = p1;			%% base options vector passed in
	i = 2;					%% start processing parameters with 2nd one
else					%% even number of parameters
	options = [				%% use defaults for base options vector
	
		%% power flow options
		1;		%% 1  - PF_ALG
		1e-8;	%% 2  - PF_TOL
		10;		%% 3  - PF_MAX_IT
		30;		%% 4  - PF_MAX_IT_FD
		0;		%% 5  - RESERVED5
		0;		%% 6  - RESERVED6
		0;		%% 7  - RESERVED7
		0;		%% 8  - RESERVED8
		0;		%% 9  - RESERVED9
		0;		%% 10 - RESERVED10
		
		%% OPF options
		0;		%% 11 - OPF_ALG_POLY
		100;	%% 12 - OPF_ALG_POLY
		200;	%% 13 - OPF_ALG_PWL
		10;		%% 14 - OPF_POLY2PWL_PTS
		0;		%% 15 - OPF_NEQ
		5e-6;	%% 16 - OPF_VIOLATION
		1e-4;	%% 17 - CONSTR_TOL_X
		1e-4;	%% 18 - CONSTR_TOL_F
		0;		%% 19 - CONSTR_MAX_IT
		3e-3;	%% 20 - LPC_TOL_GRAD
		5e-3;	%% 21 - LPC_TOL_X
		1000;	%% 22 - LPC_MAX_IT
		5;		%% 23 - LPC_MAX_RESTART
		0;		%% 24 - RESERVED24
		0;		%% 25 - RESERVED25
		0;		%% 26 - RESERVED26
		0;		%% 27 - RESERVED27
		0;		%% 28 - RESERVED28
		0;		%% 29 - RESERVED29
		0;		%% 30 - RESERVED30
		
		%% output options
		1;		%% 31 - VERBOSE
		-1;		%% 32 - OUT_ALL
		1;		%% 33 - OUT_SYS_SUM
		0;		%% 34 - OUT_AREA_SUM
		1;		%% 35 - OUT_BUS
		1;		%% 36 - OUT_BRANCH
		0;		%% 37 - OUT_GEN
		-1;		%% 38 - OUT_ALL_LIM
		1;		%% 39 - OUT_V_LIM
		1;		%% 40 - OUT_LINE_LIM
		1;		%% 41 - OUT_PG_LIM
		1;		%% 42 - OUT_QG_LIM
		0;		%% 43 - OUT_RAW
		0;		%% 44 - RESERVED44
		0;		%% 45 - RESERVED45
		0;		%% 46 - RESERVED46
		0;		%% 47 - RESERVED47
		0;		%% 48 - RESERVED48
		0;		%% 49 - RESERVED49
		0;		%% 50 - RESERVED50
	];
end

%%-----  set up option names  -----
%% power flow options
names = str2mat(	'PF_ALG', ...			%% 1
					'PF_TOL', ...			%% 2
					'PF_MAX_IT', ...		%% 3
					'PF_MAX_IT_FD', ...		%% 4
					'RESERVED5', ...		%% 5
					'RESERVED6', ...		%% 6
					'RESERVED7', ...		%% 7
					'RESERVED8', ...		%% 8
					'RESERVED9', ...		%% 9
					'RESERVED10'	);		%% 10

%% OPF options
names = str2mat(	names, ...
					'OPF_ALG', ...			%% 11
					'OPF_ALG_POLY', ...		%% 12
					'OPF_ALG_PWL', ...		%% 13
					'OPF_POLY2PWL_PTS', ...	%% 14
					'OPF_NEQ', ...			%% 15
					'OPF_VIOLATION', ...	%% 16
					'CONSTR_TOL_X', ...		%% 17
					'CONSTR_TOL_F', ...		%% 18
					'CONSTR_MAX_IT', ...	%% 19
					'LPC_TOL_GRAD'	);		%% 20
names = str2mat(	names, ...
					'LPC_TOL_X', ...		%% 21
					'LPC_MAX_IT', ...		%% 22
					'LPC_MAX_RESTART', ...	%% 23
					'RESERVED24', ...		%% 24
					'RESERVED25', ...		%% 25
					'RESERVED26', ...		%% 26
					'RESERVED27', ...		%% 27
					'RESERVED28', ...		%% 28
					'RESERVED29', ...		%% 29
					'RESERVED30'	);		%% 30

%% output options
names = str2mat(	names, ...
					'VERBOSE', ...			%% 31
					'OUT_ALL', ...			%% 32
					'OUT_SYS_SUM', ...		%% 33
					'OUT_AREA_SUM', ...		%% 34
					'OUT_BUS', ...			%% 35
					'OUT_BRANCH', ...		%% 36
					'OUT_GEN', ...			%% 37
					'OUT_ALL_LIM', ...		%% 38
					'OUT_V_LIM', ...		%% 39
					'OUT_LINE_LIM'	);		%% 40
names = str2mat(	names, ...
					'OUT_PG_LIM', ...		%% 41
					'OUT_QG_LIM', ...		%% 42
					'OUT_RAW', ...			%% 43
					'RESERVED44', ...		%% 44
					'RESERVED45', ...		%% 45
					'RESERVED46', ...		%% 46
					'RESERVED47', ...		%% 47
					'RESERVED48', ...		%% 48
					'RESERVED49', ...		%% 49
					'RESERVED50'	);		%% 50

%%-----  process parameters  -----
while i <= nargin
	%% get parameter name and value
	pname = eval(['p' int2str(i)]);
	pval = eval(['p' int2str(i+1)]);
	
	%% get parameter index
	namestr = names';
	namestr = namestr(:)';
	namelen = size(names, 2);
	pidx = ceil(findstr([pname blanks(namelen-length(pname))], namestr) / namelen);
	% fprintf('''%s'' (%d) = %d\n', pname, pidx, pval);

	%% update option
	options(pidx) = pval;

	i = i + 2;								%% go to next parameter
end

return;
