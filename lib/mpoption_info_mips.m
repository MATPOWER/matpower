function opt = mpoption_info_mips(selector)
%MPOPTION_INFO_MIPS  Returns MATPOWER option info for MIPS (optional fields).
%
%   DEFAULT_OPTS = MPOPTION_INFO_MIPS('D')
%   VALID_OPTS   = MPOPTION_INFO_MIPS('V')
%   EXCEPTIONS   = MPOPTION_INFO_MIPS('E')
%
%   Returns a structure for MIPS options for MATPOWER containing ...
%   (1) default options,
%   (2) valid options, or
%   (3) NESTED_STRUCT_COPY exceptions for setting options
%   ... depending on the value of the input argument.
%
%   This function is used by MPOPTION to set default options, check validity
%   of option names or modify option setting/copying behavior for this
%   subset of optional MATPOWER options.
%
%   See also MPOPTION.

%   MATPOWER
%   Copyright (c) 2014-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    selector = 'D';
end
switch upper(selector)
    case 'D'            %% default options
        opt = struct([]);   %% no default options
    case 'V'            %% valid options
        opt = struct(...
            'mips', struct(...
                'xi', 0.99995, ...
                'sigma', 0.1, ...
                'z0', 1, ...
                'alpha_min', 1e-8, ...
                'rho_min', 0.95, ...
                'rho_max', 1.05, ...
                'mu_threshold', 1e-5, ...
                'max_stepsize', 1e10 ...
            ) ...
        );
    case 'E'            %% exceptions used by nested_struct_copy() for applying
        opt = struct([]);   %% no exceptions
    otherwise
        error('mpoption_info_mips: ''%s'' is not a valid input argument', selector);
end
