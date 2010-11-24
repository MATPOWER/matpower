function opt = mosek_options(overrides, mpopt)
%MOSEK_OPTIONS  Sets options for MOSEK.
%
%   OPT = MOSEK_OPTIONS
%   OPT = MOSEK_OPTIONS(OVERRIDES)
%   OPT = MOSEK_OPTIONS(OVERRIDES, FNAME)
%   OPT = MOSEK_OPTIONS(OVERRIDES, MPOPT)
%
%   Sets the values for the param struct normally passed to MOSEKOPT.
%
%   Inputs are all optional, second argument must be either a string
%   (FNAME) or a vector (MPOPT):
%
%       OVERRIDES - struct containing values to override the defaults
%       FNAME - name of user-supplied function called after default
%           options are set to modify them. Calling syntax is:
%               MODIFIED_OPT = FNAME(DEFAULT_OPT);
%       MPOPT - MATPOWER options vector, uses the following entries:
%           OPF_VIOLATION (16)  - used to set opt.MSK_DPAR_INTPNT_TOL_PFEAS
%           VERBOSE (31)        - not currently used here
%           MOSEK_MAX_IT (112)  - used to set opt.MSK_IPAR_INTPNT_MAX_ITERATIONS
%           MOSEK_GAP_TOL (113) - used to set opt.MSK_DPAR_INTPNT_TOL_REL_GAP
%           MOSEK_MAX_TIME (114) - used to set opt.MSK_DPAR_OPTIMIZER_MAX_TIME
%           MOSEK_NUM_THREADS (115) - used to set opt.MSK_IPAR_INTPNT_NUM_THREADS
%           MOSEK_OPT (116)     - user option file, if MPOPT(116) is non-zero
%               non-zero it is appended to 'mosek_user_options_' to form
%               the name of a user-supplied function used as FNAME
%               described above, except with calling syntax:
%                   MODIFIED_OPT = FNAME(DEFAULT_OPT, MPOPT);
%
%   Output is a param struct to pass to MOSEKOPT.
%
%   Example:
%
%   If MPOPT(116) = 3, then after setting the default MOSEK options,
%   MOSEK_OPTIONS will execute the following user-defined function
%   to allow option overrides:
%
%       opt = mosek_user_options_3(opt, mpopt);
%
%   The contents of mosek_user_options_3.m, could be something like:
%
%       function opt = mosek_user_options_3(opt, mpopt)
%       opt.MSK_DPAR_INTPNT_TOL_DFEAS   = 1e-9;
%       opt.MSK_IPAR_SIM_MAX_ITERATIONS = 5000000;
%
%   See the Parameters reference in Appendix E of "The MOSEK
%   optimization toolbox for MATLAB manaul" for
%   details on the available options.
%
%       http://www.mosek.com/documentation/
%
%   See also MOSEKOPT, MPOPTION.

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
verbose = 2;
gaptol  = 0;
fname   = '';

%% get symbolic constant names
[r, res] = mosekopt('symbcon echo(0)');
sc = res.symbcon;

%% second argument
if nargin > 1 && ~isempty(mpopt)
    if ischar(mpopt)        %% 2nd arg is FNAME (string)
        fname = mpopt;
        have_mpopt = 0;
    else                    %% 2nd arg is MPOPT (MATPOWER options vector)
        have_mpopt = 1;
        %% (make default OPF_VIOLATION correspond to default MOSEK intpnt_tol_pfeas)
        verbose = mpopt(31);    %% VERBOSE
        if mpopt(116)           %% MOSEK_OPT
            fname = sprintf('mosek_user_options_%d', mpopt(116));
        end
    end
else
    have_mpopt = 0;
end

%%-----  set default options for MOSEK  -----
%% solution algorithm
if have_mpopt
    alg = mpopt(111);   %% MOSEK_LP_ALG
    switch alg
        case {  sc.MSK_OPTIMIZER_FREE,                  %% 0
                sc.MSK_OPTIMIZER_INTPNT,                %% 1
                sc.MSK_OPTIMIZER_PRIMAL_SIMPLEX,        %% 4
                sc.MSK_OPTIMIZER_DUAL_SIMPLEX,          %% 5
                sc.MSK_OPTIMIZER_PRIMAL_DUAL_SIMPLEX,   %% 6
                sc.MSK_OPTIMIZER_FREE_SIMPLEX,          %% 7
                sc.MSK_OPTIMIZER_CONCURRENT }           %% 10
            opt.MSK_IPAR_OPTIMIZER = alg;
        otherwise
            opt.MSK_IPAR_OPTIMIZER = sc.MSK_OPTIMIZER_FREE;
    end
    %% (make default OPF_VIOLATION correspond to default MSK_DPAR_INTPNT_TOL_PFEAS)
    opt.MSK_DPAR_INTPNT_TOL_PFEAS = mpopt(16)/500;  %% OPF_VIOLATION
    if mpopt(112)       %% MOSEK_MAX_IT
        opt.MSK_IPAR_INTPNT_MAX_ITERATIONS = mpopt(112);
    end
    if mpopt(113)       %% MOSEK_GAP_TOL
        opt.MSK_DPAR_INTPNT_TOL_REL_GAP = mpopt(113);
    end
    if mpopt(114)       %% MOSEK_MAX_TIME
        opt.MSK_DPAR_OPTIMIZER_MAX_TIME = mpopt(114);
    end
    if mpopt(115)       %% MOSEK_NUM_THREADS
        opt.MSK_IPAR_INTPNT_NUM_THREADS = mpopt(115);
    end
else
    opt.MSK_IPAR_OPTIMIZER = sc.MSK_OPTIMIZER_FREE;
end
% opt.MSK_DPAR_INTPNT_TOL_PFEAS = 1e-8;       %% primal feasibility tol
% opt.MSK_DPAR_INTPNT_TOL_DFEAS = 1e-8;       %% dual feasibility tol
% opt.MSK_DPAR_INTPNT_TOL_MU_RED = 1e-16;     %% relative complementarity gap tol
% opt.MSK_DPAR_INTPNT_TOL_REL_GAP = 1e-8;     %% relative gap termination tol
% opt.MSK_IPAR_INTPNT_MAX_ITERATIONS = 400;   %% max iterations for int point
% opt.MSK_IPAR_SIM_MAX_ITERATIONS = 10000000; %% max iterations for simplex
% opt.MSK_DPAR_OPTIMIZER_MAX_TIME = -1;       %% max time allowed (< 0 --> Inf)
% opt.MSK_IPAR_INTPNT_NUM_THREADS = 1;        %% number of threads
% opt.MSK_IPAR_PRESOLVE_USE = sc.MSK_PRESOLVE_MODE_OFF;

% if verbose == 0
%     opt.MSK_IPAR_LOG = 0;
% end

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
        opt.(names{k}) = overrides.(names{k});
    end
end
