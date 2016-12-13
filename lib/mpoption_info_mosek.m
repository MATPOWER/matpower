function opt = mpoption_info_mosek(selector)
%MPOPTION_INFO_MOSEK  Returns MATPOWER option info for MOSEK.
%
%   DEFAULT_OPTS = MPOPTION_INFO_MOSEK('D')
%   VALID_OPTS   = MPOPTION_INFO_MOSEK('V')
%   EXCEPTIONS   = MPOPTION_INFO_MOSEK('E')
%
%   Returns a structure for MOSEK options for MATPOWER containing ...
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
if have_fcn('mosek')
    switch upper(selector)
        case {'D', 'V'}     %% default and valid options
            opt = struct(...
                'mosek',        struct(...
                    'lp_alg',       0, ...
                    'max_it',       0, ...
                    'gap_tol',      0, ...
                    'max_time',     0, ...
                    'num_threads',  0, ...
                    'opts',         [], ...
                    'opt_fname',    '', ...
                    'opt',          0 ...
                ) ...
            );
        case 'E'            %% exceptions used by nested_struct_copy() for applying
            opt = struct(...
                'name',         { 'mosek.opts' }, ...
                'check',        0 ...
                );
%                 'copy_mode',    { @mosek_options } ...
        otherwise
            error('mpoption_info_mosek: ''%s'' is not a valid input argument', selector);
    end
else
    opt = struct([]);       %% MOSEK is not available
end
