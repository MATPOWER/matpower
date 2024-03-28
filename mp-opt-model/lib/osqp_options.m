function opt = osqp_options(overrides, mpopt)
% osqp_options - Sets options for OSQP.
% ::
%
%   OPT = OSQP_OPTIONS
%   OPT = OSQP_OPTIONS(OVERRIDES)
%   OPT = OSQP_OPTIONS(OVERRIDES, FNAME)
%   OPT = OSQP_OPTIONS(OVERRIDES, MPOPT)
%
%   Sets the values for the settings struct normally passed to OSQP.
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
%           osqp.opts      - struct containing values to use as OVERRIDES
%           osqp.opt_fname - name of user-supplied function used as FNAME,
%               except with calling syntax:
%                   MODIFIED_OPT = FNAME(DEFAULT_OPT, MPOPT);
%
%   Output is a settings struct to pass to OSQP.
%
%   There are multiple ways of providing values to override the default
%   options. Their precedence and order of application are as follows:
%
%   With inputs OVERRIDES and FNAME
%       1. FNAME is called
%       2. OVERRIDES are applied
%   With inputs OVERRIDES and MPOPT
%       1. FNAME (from osqp.opt_fname) is called
%       2. osqp.opts (if not empty) are applied
%       3. OVERRIDES are applied
%
%   Example:
%
%   If osqp.opt_fname = 'osqp_user_options_3', then after setting the
%   default OSQP options, OSQP_OPTIONS will execute the following
%   user-defined function to allow option overrides:
%
%       opt = osqp_user_options_3(opt, mpopt);
%
%   The contents of osqp_user_options_3.m, could be something like:
%
%       function opt = osqp_user_options_3(opt, mpopt)
%       opt.polish = 1;
%       opt.alpha  = 1;
%       opt.max_iter = 5000;
%
%   See the OSQP Solver settings page.
%
%       https://osqp.org/docs/interfaces/solver_settings.html
%
% See also osqp, mpoption.

%   MP-Opt-Model
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

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
        if isfield(mpopt.osqp, 'opt_fname') && ~isempty(mpopt.osqp.opt_fname)
            fname = mpopt.osqp.opt_fname;
        end
    end
else
    have_mpopt = 0;
end

%%-----  set default options for OSQP  -----
opt = struct( ...
    );
%         'polish', 1 ...
%         'alpha', 1 ...              %% default = 1.6 (unless problem dependent)
%         'eps_abs', 1e-9, ...        %% default = 1e-3
%         'eps_rel', 1e-9, ...        %% default = 1e-3
%         'eps_prim_inf', 1e-7, ...   %% default = 1e-4
%         'eps_dual_inf', 1e-7, ...   %% default = 1e-4
if have_mpopt
    opt.polish = 1;
    opt.eps_prim_inf = mpopt.opf.violation;
%     opt.eps_abs = 1e-4;
%     opt.eps_rel = 1e-4;
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
if have_mpopt && isfield(mpopt.osqp, 'opts') && ~isempty(mpopt.osqp.opts)
    opt = nested_struct_copy(opt, mpopt.osqp.opts);
end
if nargin > 0 && ~isempty(overrides)
    opt = nested_struct_copy(opt, overrides);
end
