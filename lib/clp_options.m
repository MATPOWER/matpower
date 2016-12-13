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
%   Copyright (c) 2010-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

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
