function opt = mpoption_info_clp(selector)
%MPOPTION_INFO_CLP  Returns MATPOWER option info for CLP.
%
%   DEFAULT_OPTS = MPOPTION_INFO_CLP('D')
%   VALID_OPTS   = MPOPTION_INFO_CLP('V')
%   EXCEPTIONS   = MPOPTION_INFO_CLP('E')
%
%   Returns a structure for CLP options for MATPOWER containing ...
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
%   Copyright (c) 2014-2015, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    selector = 'D';
end
if have_fcn('clp')
    switch upper(selector)
        case {'D', 'V'}     %% default and valid options
            opt = struct(...
                'clp',           struct(...
                    'opts',            [], ...
                    'opt_fname',       '' ...
                ) ...
            );
        case 'E'            %% exceptions used by nested_struct_copy() for applying
            opt = struct(...
                'name',         { 'clp.opts' }, ...
                'check',        0 ...
                );
%                 'copy_mode',    { @clp_options } ...
        otherwise
            error('mpoption_info_clp: ''%s'' is not a valid input argument', selector);
    end
else
    opt = struct([]);       %% CLP is not available
end
