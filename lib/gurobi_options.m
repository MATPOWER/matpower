function opt = gurobi_options(overrides, mpopt)
%GUROBI_OPTIONS  Sets options for GUROBI.
%
%   OPT = GUROBI_OPTIONS
%   OPT = GUROBI_OPTIONS(OVERRIDES)
%   OPT = GUROBI_OPTIONS(OVERRIDES, FNAME)
%   OPT = GUROBI_OPTIONS(OVERRIDES, MPOPT)
%
%   Sets the values for the options struct normally passed to GUROBI_MEX.
%
%   Inputs are all optional, second argument must be either a string
%   (FNAME) or a vector (MPOPT):
%
%       OVERRIDES - struct containing values to override the defaults
%       FNAME - name of user-supplied function called after default
%           options are set to modify them. Calling syntax is:
%               MODIFIED_OPT = FNAME(DEFAULT_OPT);
%       MPOPT - MATPOWER options vector, uses the following entries:
%           OPF_VIOLATION (16)  - used to set opt.FeasibilityTol
%           VERBOSE (31)        - used to set opt.DisplayInterval, opt.Display
%           GRB_METHOD (121)    - used to set opt.Method
%           GRB_TIMELIMIT (122) - used to set opt.TimeLimit (seconds)
%           GRB_THREADS (123)   - used to set opt.Threads
%           GRB_OPT (124)       - user option file, if MPOPT(124) is non-zero
%               it is appended to 'gurobi_user_options_' to form the name of a
%               user-supplied function used as FNAME described above, except
%               with calling syntax:
%                   MODIFIED_OPT = FNAME(DEFAULT_OPT, MPOPT);
%
%   Output is an options struct to pass to GUROBI_MEX.
%
%   Example:
%
%   If MPOPT(124) = 3, then after setting the default GUROBI options,
%   GUROBI_OPTIONS will execute the following user-defined function
%   to allow option overrides:
%
%       opt = gurobi_user_options_3(opt, mpopt);
%
%   The contents of gurobi_user_options_3.m, could be something like:
%
%       function opt = gurobi_user_options_3(opt, mpopt)
%       opt.OptimalityTol   = 1e-9;
%       opt.IterationLimit  = 3000;
%       opt.BarIterLimit    = 200;
%       opt.Crossover       = 0;
%       opt.Presolve        = 0;
%
%   For details on the available options, see the "Parameters" section
%   of the "Gurobi Optimizer Reference Manual" at:
%
%       http://www.gurobi.com/doc/45/refman/
%
%   See also GUROBI_MEX, MPOPTION.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2010-2011 by Power System Engineering Research Center (PSERC)
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
fname   = '';

%% second argument
if nargin > 1 && ~isempty(mpopt)
    if ischar(mpopt)        %% 2nd arg is FNAME (string)
        fname = mpopt;
        have_mpopt = 0;
    else                    %% 2nd arg is MPOPT (MATPOWER options vector)
        have_mpopt = 1;
        verbose = mpopt(31);    %% VERBOSE
        if mpopt(124)           %% GRB_OPT
            fname = sprintf('gurobi_user_options_%d', mpopt(124));
        end
    end
else
    have_mpopt = 0;
end

%%-----  set default options for CPLEX  -----
% opt.OptimalityTol = 1e-6;
% opt.Presolve = -1;              %% -1 - auto, 0 - no, 1 - conserv, 2 - aggressive=
% opt.LogFile = 'qps_gurobi.log';
if have_mpopt
    %% (make default OPF_VIOLATION correspond to default FeasibilityTol)
    opt.FeasibilityTol  = mpopt(16)/5;  %% OPF_VIOLATION
    opt.Method          = mpopt(121);   %% GRB_METHOD
    opt.TimeLimit       = mpopt(122);   %% GRB_TIMELIMIT
    opt.Threads         = mpopt(123);   %% GRB_THREADS
else
    opt.Method          = 1;            %% dual simplex
end
opt.Display = min(verbose, 3);
if verbose
    opt.DisplayInterval = 1;
else
    opt.DisplayInterval = Inf;
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
        opt.(names{k}) = overrides.(names{k});
    end
end
