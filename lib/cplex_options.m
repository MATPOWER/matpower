function opt = cplex_options(overrides, mpopt)
%CPLEX_OPTIONS  Sets options for CPLEX.
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
%   (FNAME) or a vector (MPOPT):
%
%       OVERRIDES - struct containing values to override the defaults
%       FNAME - name of user-supplied function called after default
%           options are set to modify them. Calling syntax is:
%               MODIFIED_OPT = FNAME(DEFAULT_OPT);
%       MPOPT - MATPOWER options vector, used to set:
%           feasibility tol  - based on MPOPT(16) (OPF_VIOLATION)
%           output verbosity - based on MPOPT(31) (VERBOSE)
%           LP solver alg    - based on MPOPT(95) (CPLEX_LPMETHOD)
%           QP solver alg    - based on MPOPT(96) (CPLEX_QPMETHOD)
%           user option file - based on MPOPT(97) (CPLEX_OPT)
%               If MPOPT(97) (CPLEX_OPT) is non-zero it is appended to
%               'cplex_user_options_' to form the name of a
%               user-supplied function used as FNAME described above,
%               except with calling syntax:
%               MODIFIED_OPT = FNAME(DEFAULT_OPT, MPOPT);
%
%   Output is an options struct to pass to CPLEXOPTIMSET.
%
%   Example:
%
%   If MPOPT(97) = 3, then after setting the default CPLEX options,
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
%   For details on the available options, see the "Parameters Reference
%   Manual" section of the CPLEX documentation at:
%
%       http://publib.boulder.ibm.com/infocenter/cosinfoc/v12r2/
%
%   See also CPLEXLP, CPLEXQP, MPOPTION.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2010 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

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
    else                    %% 2nd arg is MPOPT (MATPOWER options vector)
        have_mpopt = 1;
        %% (make default OPF_VIOLATION correspond to default CPLEX feastol)
        feastol = mpopt(16)/5;  %% OPF_VIOLATION
        verbose = mpopt(31);    %% VERBOSE
        lpmethod = mpopt(95);   %% CPLEX_LPMETHOD
        qpmethod = mpopt(96);   %% CPLEX_QPMETHOD
        if mpopt(97)            %% CPLEX_OPT
            fname = sprintf('cplex_user_options_%d', mpopt(97));
        end
    end
else
    have_mpopt = 0;
end

%%-----  set default options for CPLEX  -----
opt = cplexoptimset('cplex');
opt.simplex.tolerances.feasibility = feastol;

%% printing
vrb = max([0 verbose-1]);
cplex_opt.barrier.display   = vrb;
cplex_opt.conflict.display  = vrb;
cplex_opt.mip.display       = vrb;
cplex_opt.sifting.display   = vrb;
cplex_opt.simplex.display   = vrb;
cplex_opt.tune.display      = vrb;

%% solution algorithm
if have_mpopt
    opt.lpmethod = mpopt(95);   %% CPLEX_LPMETHOD
    opt.qpmethod = mpopt(96);   %% CPLEX_QPMETHOD
% else
%     opt.lpmethod = 2;
%     opt.qpmethod = 2;
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
if nargin > 0 && ~isempty(overrides)
    names = fieldnames(overrides);
    for k = 1:length(names)
        if isstruct(overrides.(names{k}))
            names2 = fieldnames(overrides.(names{k}));
            for k2 = 1:length(names2)
                if isstruct(overrides.(names{k}).(names2{k2}))
                    names3 = fieldnames(overrides.(names{k}).(names2{k2}));
                    for k3 = 1:length(names3)
                        opt.(names{k}).(names2{k2}).(names3{k3}) = overrides.(names{k}).(names2{k2}).(names3{k3});
                    end
                else
                    opt.(names{k}).(names2{k2}) = overrides.(names{k}).(names2{k2});
                end
            end
        else
            opt.(names{k}) = overrides.(names{k});
        end
    end
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
%               iteration:      2.1000e+09
%               objrange:       1.0000e+20
%           ordering:       0
%           qcpconvergetol: 1.0000e-07
%           startalg:       1
%       clocktype:      2
%       conflict:       [1x1 struct]
%           display:        1
%       emphasis:       [1x1 struct]
%           memory:         0
%           mip:            0
%           numerical:      0
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
%               nodes:          2.1000e+09
%               polishtime:     0
%               populate:       20
%               probetime:      1.0000e+75
%               repairtries:    0
%               solutions:      2.1000e+09
%               strongcand:     10
%               strongit:       0
%               submipnodelim:  500
%               treememory:     1.0000e+75
%           ordertype:      0
%           polishafter:    [1x1 struct]
%               absmipgap:      0
%               mipgap:         0
%               nodes:          2.1000e+09
%               solutions:      2.1000e+09
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
%       parallel:       0
%       preprocessing:  [1x1 struct]
%           aggregator:     -1
%           boundstrength:  -1
%           coeffreduce:    2
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
%           constraints:    30000
%           datacheck:      0
%           nonzeros:       250000
%           qpnonzeros:     5000
%           scale:          0
%           variables:      60000
%       sifting:        [1x1 struct]
%           algorithm:      0
%           display:        1
%           iterations:     2.1000e+09
%       simplex:        [1x1 struct]
%           crash:          1
%           dgradient:      0
%           display:        1
%           limits:         [1x1 struct]
%               iterations:     2.1000e+09
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
%       threads:        0
%       timelimit:      1.0000e+75
%       tune:           [1x1 struct]
%           display:        1
%           measure:        1
%           repeat:         1
%           timelimit:      10000
%       workdir:        '.'
%       workmem:        128
