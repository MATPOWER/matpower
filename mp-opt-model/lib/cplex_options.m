function opt = cplex_options(overrides, mpopt)
% cplex_options - Sets options for CPLEX.
% ::
%
%   OPT = CPLEX_OPTIONS
%   OPT = CPLEX_OPTIONS(OVERRIDES)
%   OPT = CPLEX_OPTIONS(OVERRIDES, FNAME)
%   OPT = CPLEX_OPTIONS(OVERRIDES, MPOPT)
%
%   Sets the values for the options struct normally passed to
%   CPLEXOPTIMSET.
%
%   Inputs are all optional, second argument must be either a string
%   (FNAME) or a struct (MPOPT):
%
%       OVERRIDES - struct containing values to override the defaults
%       FNAME - name of user-supplied function called after default
%           options are set to modify them. Calling syntax is:
%                   MODIFIED_OPT = FNAME(DEFAULT_OPT);
%       MPOPT - MATPOWER options struct, uses the following fields:
%           opf.violation    - used to set opt.simplex.tolerances.feasibility
%           verbose          - used to set opt.barrier.display,
%               opt.conflict.display, opt.mip.display, opt.sifting.display,
%               opt.simplex.display, opt.tune.display
%           cplex.lpmethod   - used to set opt.lpmethod
%           cplex.qpmethod   - used to set opt.qpmethod
%           cplex.opts       - struct containing values to use as OVERRIDES
%           cplex.opt_fname  - name of user-supplied function used as FNAME,
%               except with calling syntax:
%                   MODIFIED_OPT = FNAME(DEFAULT_OPT, MPOPT);
%           cplex.opt        - numbered user option function, if and only if
%               cplex.opt_fname is empty and cplex.opt is non-zero, the value
%               of cplex.opt_fname is generated by appending cplex.opt to
%               'cplex_user_options_' (for backward compatibility with old
%               MATPOWER option CPLEX_OPT).
%
%   Output is an options struct to pass to CPLEXOPTIMSET.
%
%   There are multiple ways of providing values to override the default
%   options. Their precedence and order of application are as follows:
%
%   With inputs OVERRIDES and FNAME
%       1. FNAME is called
%       2. OVERRIDES are applied
%   With inputs OVERRIDES and MPOPT
%       1. FNAME (from cplex.opt_fname or cplex.opt) is called
%       2. cplex.opts (if not empty) are applied
%       3. OVERRIDES are applied
%
%   Example:
%
%   If cplex.opt = 3, then after setting the default CPLEX options,
%   CPLEX_OPTIONS will execute the following user-defined function
%   to allow option overrides:
%
%       opt = cplex_user_options_3(opt, mpopt);
%
%   The contents of cplex_user_options_3.m, could be something like:
%
%       function opt = cplex_user_options_3(opt, mpopt)
%       opt.threads          = 2;
%       opt.simplex.refactor = 1;
%       opt.timelimit        = 10000;
%
%   For details on the available options, see the "Parameters of CPLEX"
%   section of the CPLEX documentation at:
%
%       http://www.ibm.com/support/knowledgecenter/SSSA5P
%
% See also cplexlp, cplexqp, mpoption.

%   MP-Opt-Model
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%%-----  initialization and arg handling  -----
%% defaults
verbose = 1;
feastol = 1e-6;
fname   = '';

%% second argument
if nargin > 1 && ~isempty(mpopt)
    if ischar(mpopt)        %% 2nd arg is FNAME (string)
        fname = mpopt;
        have_mpopt = 0;
    else                    %% 2nd arg is MPOPT (MATPOWER options struct)
        have_mpopt = 1;
        %% (make default opf.violation correspond to default CPLEX feastol)
        feastol  = mpopt.opf.violation/5;
        verbose  = mpopt.verbose;
        if isfield(mpopt.cplex, 'opt_fname') && ~isempty(mpopt.cplex.opt_fname)
            fname = mpopt.cplex.opt_fname;
        elseif mpopt.cplex.opt
            fname = sprintf('cplex_user_options_%d', mpopt.cplex.opt);
        end
    end
else
    have_mpopt = 0;
end

%%-----  set default options for CPLEX  -----
opt = struct( ...
    'simplex', struct('tolerances', struct('feasibility', feastol)), ...
    'output', struct('clonelog', -1) ...
);

%% printing
if verbose > 2
    opt.display = 'iter';
elseif verbose > 1
    opt.display = 'on';
elseif verbose > 0
    opt.display = 'off';
end

%% solution algorithm
if have_mpopt
    opt.lpmethod = mpopt.cplex.lpmethod;
    opt.qpmethod = mpopt.cplex.qpmethod;
end

%%-----  call user function to modify defaults  -----
if ~isempty(fname)
    if have_mpopt
        opt = feval(fname, opt, mpopt);
    else
        opt = feval(fname, opt);
    end
end

%%-----  apply overrides  -----
if have_mpopt && isfield(mpopt.cplex, 'opts') && ~isempty(mpopt.cplex.opts)
    opt = nested_struct_copy(opt, mpopt.cplex.opts);
end
if nargin > 0 && ~isempty(overrides)
    opt = nested_struct_copy(opt, overrides);
end


%--------------------------  Default Options Struct  --------------------------
% as returned by ...
%   >> opt = cplexoptimset('cplex')
%
%   opt =
%       advance:        1
%       barrier:        [1x1 struct]
%           algorithm:      0
%           colnonzeros:    0
%           convergetol:    1.0000e-08
%           crossover:      0
%           display:        1
%           limits:         [1x1 struct]
%               corrections:    -1
%               growth:         1.0000e+12
%               iteration:      9.2234e+18
%               objrange:       1.0000e+20
%           ordering:       0
%           qcpconvergetol: 1.0000e-07
%           startalg:       1
%       clocktype:      2
%       conflict:       [1x1 struct]
%           display:        1
%       diagnostics:    'off'
%       emphasis:       [1x1 struct]
%           memory:         0
%           mip:            0
%           numerical:      0
%       exportmodel:    ''
%       feasopt:        [1x1 struct]
%           mode:           0
%           tolerance:      1.0000e-06
%       lpmethod:       0
%       mip:            [1x1 struct]
%           cuts:           [1x1 struct]
%               cliques:        0
%               covers:         0
%               disjunctive:    0
%               flowcovers:     0
%               gomory:         0
%               gubcovers:      0
%               implied:        0
%               mcfcut:         0
%               mircut:         0
%               pathcut:        0
%               zerohalfcut:    0
%           display:        2
%           interval:       0
%           limits:         [1x1 struct]
%               aggforcut:      3
%               auxrootthreads: 0
%               cutpasses:      0
%               cutsfactor:     4
%               eachcutlimit:   2.1000e+09
%               gomorycand:     200
%               gomorypass:     0
%               nodes:          9.2234e+18
%               polishtime:     0
%               populate:       20
%               probetime:      1.0000e+75
%               repairtries:    0
%               solutions:      9.2234e+18
%               strongcand:     10
%               strongit:       0
%               submipnodelim:  500
%               treememory:     1.0000e+75
%           ordertype:      0
%           polishafter:    [1x1 struct]
%               absmipgap:      0
%               mipgap:         0
%               nodes:          9.2234e+18
%               solutions:      9.2234e+18
%               time:           1.0000e+75
%           pool:           [1x1 struct]
%               absgap:         1.0000e+75
%               capacity:       2.1000e+09
%               intensity:      0
%               relgap:         1.0000e+75
%               replace:        0
%           strategy:       [1x1 struct]
%               backtrack:      0.9999
%               bbinterval:     7
%               branch:         0
%               dive:           0
%               file:           1
%               fpheur:         0
%               heuristicfreq:  0
%               kappastats:     0
%               lbheur:         0
%               miqcpstrat:     0
%               nodeselect:     1
%               order:          1
%               presolvenode:   0
%               probe:          0
%               rinsheur:       0
%               search:         0
%               startalgorithm: 0
%               subalgorithm:   0
%               variableselect: 0
%           tolerances:     [1x1 struct]
%               absmipgap:      1.0000e-06
%               integrality:    1.0000e-05
%               lowercutoff:    -1.0000e+75
%               mipgap:         1.0000e-04
%               objdifference:  0
%               relobjdifference: 0
%               uppercutoff:    1.0000e+75
%       output:         [1x1 struct]
%           clonelog:       1
%           intsolfileprefix: ''
%           mpslong:        1
%           writelevel:     0
%       parallel:       0
%       preprocessing:  [1x1 struct]
%           aggregator:     -1
%           boundstrength:  -1
%           coeffreduce:    -1
%           dependency:     -1
%           dual:           0
%           fill:           10
%           linear:         1
%           numpass:        -1
%           presolve:       1
%           qpmakepsd:      1
%           reduce:         3
%           relax:          -1
%           repeatpresolve: -1
%           symmetry:       -1
%       qpmethod:       0
%       read:           [1x1 struct]
%           apiencoding:    ''
%           constraints:    30000
%           datacheck:      0
%           fileencoding:   'ISO-8859-1'
%           nonzeros:       250000
%           qpnonzeros:     5000
%           scale:          0
%           variables:      60000
%       sifting:        [1x1 struct]
%           algorithm:      0
%           display:        1
%           iterations:     9.2234e+18
%       simplex:        [1x1 struct]
%           crash:          1
%           dgradient:      0
%           display:        1
%           limits:         [1x1 struct]
%               iterations:     9.2234e+18
%               lowerobj:       -1.0000e+75
%               perturbation:   0
%               singularity:    10
%               upperobj:       1.0000e+75
%           perturbation:   [1x1 struct]
%               indicator:      0
%               constant:       1.0000e-06
%           pgradient:      0
%           pricing:        0
%           refactor:       0
%           tolerances:     [1x1 struct]
%               feasibility:    1.0000e-06
%               markowitz:      0.0100
%               optimality:     1.0000e-06
%       solutiontarget: 0
%       threads:        0
%       timelimit:      1.0000e+75
%       tune:           [1x1 struct]
%           display:        1
%           measure:        1
%           repeat:         1
%           timelimit:      10000
%       workdir:        '.'
%       workmem:        128
