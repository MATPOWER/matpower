function opt = ipopt_options(overrides, mpopt)
%IPOPT_OPTIONS  Sets default options for IPOPT.
%
%   OPT = IPOPT_OPTIONS
%   OPT = IPOPT_OPTIONS(OVERRIDES)
%   OPT = IPOPT_OPTIONS(OVERRIDES, FNAME)
%   OPT = IPOPT_OPTIONS(OVERRIDES, MPOPT)
%
%   Sets the default values for the options.ipopt struct normally
%   passed to IPOPT.
%
%   Inputs are all optional, second argument must be either a string
%   (FNAME) or a vector (MPOPT):
%
%       OVERRIDES - struct containing values to override the defaults
%       FNAME - name of user-supplied function called after default
%           options are set to modify them. Calling syntax is:
%               MODIFIED_OPT = FNAME(DEFAULT_OPT);
%       MPOPT - MATPOWER options vector, used to set print_level
%           based on MPOPT(31) (VERBOSE). If MPOPT(60) (IPOPT_OPT)
%           is non-zero it is appended to 'ipopt_user_options_' to
%           form the name of a user-supplied function used as FNAME
%           described above, except with calling syntax:
%               MODIFIED_OPT = FNAME(DEFAULT_OPT, MPOPT);
%
%   Output is an options.ipopt struct to pass to IPOPT.
%
%   Example:
%
%   If MPOPT(60) = 3, then after setting the default IPOPT options,
%   IPOPT_OPTIONS will execute the following user-defined function
%   to allow option overrides:
%
%       opt = ipopt_user_options_3(opt, mpopt);
%
%   The contents of ipopt_user_options_3.m, could be something like:
%
%       function opt = ipopt_user_options_3(opt, mpopt)
%       opt.nlp_scaling_method = 'none';
%       opt.max_iter           = 500;
%       opt.derivative_test    = 'first-order';
%
%   See the options reference section in the IPOPT documentation for
%   details on the available options.
%
%       http://www.coin-or.org/Ipopt/documentation/
%
%   See also IPOPT, MPOPTION.

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

%% defaults
verbose = 2;
fname   = '';

%% second argument
if nargin > 1 && ~isempty(mpopt)
    if ischar(mpopt)        %% 2nd arg is FNAME (string)
        fname = mpopt;
    else                    %% 2nd arg is MPOPT (MATPOWER options vector)
        verbose = mpopt(31);    %% VERBOSE
        if mpopt(60)            %% IPOPT_OPT
            fname = sprintf('ipopt_user_options_%d', mpopt(60));
        end
    end
end

%% set default options for ipopt
if verbose
    opt.print_level = min(12, verbose*2+1);
else
    opt.print_level = 0;
end
%opt.print_options_documentation = 'yes';
%opt.nlp_scaling_method = 'none';
opt.mu_strategy  = 'adaptive';
%opt.tol          = 1e-7;
opt.tol           = 1e-12;
opt.max_iter      = 200;

%% call user function to modify defaults
if ~isempty(fname)
    if ischar(mpopt)        %% 2nd arg is FNAME (string)
        opt = feval(fname, opt);
    else                    %% 2nd arg is MPOPT (MATPOWER options vector)
        opt = feval(fname, opt, mpopt);
    end
end

%% apply overrides
if nargin > 0 && ~isempty(overrides)
    names = fieldnames(overrides);
    for k = 1:length(names)
        opt.(names{k}) = overrides.(names{k});
    end
end
