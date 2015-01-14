function opt = clp_options(overrides, mpopt)
%CLP_OPTIONS  Sets options for CLP.
%
%   OPT = CLP_OPTIONS
%   OPT = CLP_OPTIONS(OVERRIDES)
%   OPT = CLP_OPTIONS(OVERRIDES, FNAME)
%   OPT = CLP_OPTIONS(OVERRIDES, MPOPT)
%
%   Sets the values for the options struct normally passed to CLP.
%
%   Inputs are all optional, second argument must be either a string
%   (FNAME) or a struct (MPOPT):
%
%       OVERRIDES - struct containing values to override the defaults
%       FNAME - name of user-supplied function called after default
%           options are set to modify them. Calling syntax is:
%                   MODIFIED_OPT = FNAME(DEFAULT_OPT);
%       MPOPT - MATPOWER options struct, uses the following fields:
%           verbose        - used to set opt.msglev
%           clp.opts       - struct containing values to use as OVERRIDES
%           clp.opt_fname  - name of user-supplied function used as FNAME,
%               except with calling syntax:
%                   MODIFIED_OPT = FNAME(DEFAULT_OPT, MPOPT);
%
%   Output is an options struct to pass to CLP.
%
%   There are multiple ways of providing values to override the default
%   options. Their precedence and order of application are as follows:
%
%   With inputs OVERRIDES and FNAME
%       1. FNAME is called
%       2. OVERRIDES are applied
%   With inputs OVERRIDES and MPOPT
%       1. FNAME (from clp.opt_fname) is called
%       2. clp.opts (if not empty) are applied
%       3. OVERRIDES are applied
%
%   Example:
%
%   If clp.opt_fname = 'clp_user_options_3', then after setting the
%   default CLP options, CLP_OPTIONS will execute the following
%   user-defined function to allow option overrides:
%
%       opt = clp_user_options_3(opt, mpopt);
%
%   The contents of clp_user_options_3.m, could be something like:
%
%       function opt = clp_user_options_3(opt, mpopt)
%       opt.algorithm   = 1;
%       opt.numThreads  = 2;
%
%   See the documentation for the CLP Options by typing 'help clp'.
%
%   See also CLP, MPOPTION.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2010-2015 by Power System Engineering Research Center (PSERC)
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
fname   = '';

%% second argument
if nargin > 1 && ~isempty(mpopt)
    if ischar(mpopt)        %% 2nd arg is FNAME (string)
        fname = mpopt;
        have_mpopt = 0;
    else                    %% 2nd arg is MPOPT (MATPOWER options struct)
        have_mpopt = 1;
        verbose = mpopt.verbose;
        if isfield(mpopt.clp, 'opt_fname') && ~isempty(mpopt.clp.opt_fname)
            fname = mpopt.clp.opt_fname;
        end
    end
else
    have_mpopt = 0;
end

%%-----  set default options for CLP  -----
%% printing
if have_fcn('opti_clp')
    opt.display = verbose;
else
    opt.verbose = verbose;
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
if have_mpopt && isfield(mpopt.clp, 'opts') && ~isempty(mpopt.clp.opts)
    opt = nested_struct_copy(opt, mpopt.clp.opts);
end
if nargin > 0 && ~isempty(overrides)
    opt = nested_struct_copy(opt, overrides);
end
