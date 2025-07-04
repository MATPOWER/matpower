function opt = artelys_knitro_options(overrides, mpopt)
% artelys_knitro_options - Sets options for Artelys Knitro (version 13.x and greater).
% ::
%
%   OPT = ARTELYS_KNITRO_OPTIONS
%   OPT = ARTELYS_KNITRO_OPTIONS(OVERRIDES)
%   OPT = ARTELYS_KNITRO_OPTIONS(OVERRIDES, FNAME)
%   OPT = ARTELYS_KNITRO_OPTIONS(OVERRIDES, MPOPT)
%
%   Sets the values for the options struct normally passed to Artelys Knitro.
%
%   Inputs are all optional, second argument must be either a string
%   (FNAME) or a struct (MPOPT):
%
%       OVERRIDES - struct containing values to override the defaults
%       FNAME - name of user-supplied function called after default
%           options are set to modify them. Calling syntax is:
%                   MODIFIED_OPT = FNAME(DEFAULT_OPT);
%       MPOPT - MATPOWER options struct, uses the following fields:
%           opf.violation    - used to set opt.FeasibilityTol
%           verbose          - used to set opt.DisplayInterval,
%                                 opt.OutputFlag, opt.LogToConsole
%           knitro.tol_x     - used to set opt.xtol
%           knitro.tol_f     - used to set opt.ftol
%           knitro.maxit     - used to set opt.maxit
%           knitro.opts      - struct containing values to use as OVERRIDES
%           knitro.opt_fname - name of user-supplied options file, to be
%               passed to the solver as a 'knitroOptsFile' input argument,
%               overriding any options it specifies
%           knitro.opt      - numbered user option file, if and only if
%               knitro.opt_fname is empty and knitro.opt is a non-zero
%               integer N then knitro.opt_fname is auto generated as:
%               'knitro_user_options_N.txt'
%
%   Output is a parameter struct to pass to Artelys Knitro.
%
%   There are multiple ways of providing values to override the default
%   options. Their precedence and order of application are as follows:
%
%   With inputs OVERRIDES and FNAME
%       1. FNAME is called
%       2. OVERRIDES are applied
%   With inputs OVERRIDES and MPOPT
%       1. knitro.opts (if not empty) are applied
%       2. OVERRIDES are applied
%       3. Options in options file (if specified) are applied
%
%   For details on the available options, see the "Knitro user options"
%   section of the "Artelys Knitro User's Manual" at:
%
%       https://www.artelys.com/app/docs/knitro/index.html
%
% See also knitro_options, mpoption.

%   MP-Opt-Model
%   Copyright (c) 2010-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%   and Ray Zimmerman, PSERC Cornell
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
    end
else
    have_mpopt = 0;
end

%%-----  set default options for Knitro  -----
if have_feature('knitro')
    if have_feature('knitro', 'vnum') >= 13
        opt = knitro_options;
%         opt.bar_directinterval = 0; %% previously set to 0 by default, since
%                                     %% it can help some problems to solve
    else
        kv = knitrover;
        error('artelys_knitro_options: requires Aretelys Knitro version 13.x or later (installed version is %s)', kv)
    end
else
    error('artelys_knitro_options: requires Aretelys Knitro version 13.x or later');
end

if have_mpopt
    %% (make default opf.violation correspond to default feastol)
    opt.feastol_abs = mpopt.opf.violation;
    opt.xtol     = mpopt.knitro.tol_x;
    opt.ftol     = mpopt.knitro.tol_f;
    opt.maxit    = mpopt.knitro.maxit;
end
if verbose > 1
    opt.outlev = 3;
    if verbose > 2
        opt.outlev = 4;
    end
else
    opt.outlev = 0;
end

%%-----  call user function to modify defaults  -----
if ~isempty(fname)
    opt = feval(fname, opt);
end

%%-----  apply overrides  -----
if have_mpopt && isfield(mpopt.knitro, 'opts') && ~isempty(mpopt.knitro.opts)
    opt = nested_struct_copy(opt, mpopt.knitro.opts);
end
if nargin > 0 && ~isempty(overrides)
    opt = nested_struct_copy(opt, overrides);
end
