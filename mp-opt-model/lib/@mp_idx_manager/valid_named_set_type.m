function str = valid_named_set_type(obj, set_type)
% valid_named_set_type - Returns a label for the given named set type or empty.
% ::
%
%   -----  PRIVATE METHOD  -----
%
%   STR = OBJ.VALID_NAMED_SET_TYPE(SET_TYPE)
%
%   Returns a string corresponding to the type of the named set for
%   valid values of SET_TYPE, otherwise an empty string. Can be used to
%   check whether SET_TYPE is a valid type.
%
%   Valid set types, and the corresponding values are defined as a struct
%   returned by the DEF_SET_TYPES method.
%
%   E.g. The following method ...
%       function obj = def_set_types(obj)
%           obj.set_types = struct(...
%                   'var', 'variable', ...
%                   'lin', 'linear constraint' ...
%               );
%       end
%
%   ... indicates the following valid SET_TYPES and corresponding
%   return values:
%
%       SET_TYPE    STR
%       'var'       'variable'
%       'lin'       'linear constraint'

%   MP-Opt-Model
%   Copyright (c) 2017-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if isfield(obj.set_types, set_type)
    str = obj.set_types.(set_type);
else
    str = '';
end
