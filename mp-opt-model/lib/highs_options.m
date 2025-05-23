function opt = highs_options(overrides, mpopt)
% highs_options - Sets options for HiGHS.
% ::
%
%   OPT = HIGHS_OPTIONS
%   OPT = HIGHS_OPTIONS(OVERRIDES)
%   OPT = HIGHS_OPTIONS(OVERRIDES, FNAME)
%   OPT = HIGHS_OPTIONS(OVERRIDES, MPOPT)
%
%   Sets the values for the settings struct normally passed to HiGHS.
%
%   Inputs are all optional, second argument must be either a string
%   (FNAME) or a struct (MPOPT):
%
%       OVERRIDES - struct containing values to override the defaults
%       FNAME - name of user-supplied function called after default
%           options are set to modify them. Calling syntax is:
%                   MODIFIED_OPT = FNAME(DEFAULT_OPT);
%       MPOPT - MATPOWER options struct, uses the following fields:
%           opf.violation  - used to set opt.eps_prim_inf
%           verbose        - used to set opt.verbose
%           highs.opts     - struct containing values to use as OVERRIDES
%           highs.opt_fname - name of user-supplied function used as FNAME,
%               except with calling syntax:
%                   MODIFIED_OPT = FNAME(DEFAULT_OPT, MPOPT);
%
%   Output is a settings struct to pass to HiGHS.
%
%   There are multiple ways of providing values to override the default
%   options. Their precedence and order of application are as follows:
%
%   With inputs OVERRIDES and FNAME
%       1. FNAME is called
%       2. OVERRIDES are applied
%   With inputs OVERRIDES and MPOPT
%       1. FNAME (from highs.opt_fname) is called
%       2. highs.opts (if not empty) are applied
%       3. OVERRIDES are applied
%
%   Example:
%
%   If highs.opt_fname = 'highs_user_options_3', then after setting the
%   default HiGHS options, HIGHS_OPTIONS will execute the following
%   user-defined function to allow option overrides:
%
%       opt = highs_user_options_3(opt, mpopt);
%
%   The contents of highs_user_options_3.m, could be something like:
%
%       function opt = highs_user_options_3(opt, mpopt)
%       opt.solver = 'ipm';
%       opt.dual_feasibility_tolerance  = 1e-10;
%       opt.primal_feasibility_tolerance  = 1e-10;
%       opt.threads = 2;
%
%   See the HiGHS Solver settings page.
%
%       https://ergo-code.github.io/HiGHS/dev/options/definitions/
%
% See also callhighs, mpoption.

%   MP-Opt-Model
%   Copyright (c) 2010-2025, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%%-----  initialization and arg handling  -----
%% defaults
verbose = 1;
fname   = '';

%% second argument
if nargin > 1 && ~isempty(mpopt)
    if ischar(mpopt)        %% 2nd arg is FNAME (string)
        fname = mpopt;
        have_mpopt = 0;
    else                    %% 2nd arg is MPOPT (MATPOWER options struct)
        have_mpopt = 1;
        verbose = mpopt.verbose;
        if isfield(mpopt.highs, 'opt_fname') && ~isempty(mpopt.highs.opt_fname)
            fname = mpopt.highs.opt_fname;
        end
    end
else
    have_mpopt = 0;
end

%%-----  set default options for HiGHS  -----
opt = struct();
if have_mpopt
    opt.primal_feasibility_tolerance = mpopt.opf.violation;
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
if have_mpopt && isfield(mpopt.highs, 'opts') && ~isempty(mpopt.highs.opts)
    opt = nested_struct_copy(opt, mpopt.highs.opts);
end
if nargin > 0 && ~isempty(overrides)
    opt = nested_struct_copy(opt, overrides);
end
