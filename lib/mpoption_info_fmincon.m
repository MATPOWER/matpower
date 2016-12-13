function opt = mpoption_info_fmincon(selector)
%MPOPTION_INFO_FMINCON  Returns MATPOWER option info for FMINCON.
%
%   DEFAULT_OPTS = MPOPTION_INFO_FMINCON('D')
%   VALID_OPTS   = MPOPTION_INFO_FMINCON('V')
%   EXCEPTIONS   = MPOPTION_INFO_FMINCON('E')
%
%   Returns a structure for FMINCON options for MATPOWER containing ...
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
if have_fcn('fmincon')
    switch upper(selector)
        case {'D', 'V'}     %% default and valid options
            opt = struct(...
                'fmincon',  struct(...
                    'alg',      4, ...
                    'tol_x',    1e-4, ...
                    'tol_f',    1e-4, ...
                    'max_it',   0 ...
                ) ...  %         'opt_fname', '', 'opts', []
            );
        case 'E'            %% exceptions used by nested_struct_copy() for applying
            opt = struct([]);   %% no exceptions
%             opt = struct(...
%                 'name',         { 'fmincon.opts' }, ...
%                 'check',        0, ...
%                 'copy_mode',    { '' } ...
%                 );
        otherwise
            error('mpoption_info_fmincon: ''%s'' is not a valid input argument', selector);
    end
else
    opt = struct([]);       %% FMINCON is not available
end
