function val = get(obj, varargin)
%GET  Returns the value of a field.
%   VAL = OBJ.GET(FIELD1, FIELD2, ...)
%
%   Example:
%       var_order = obj.get('var', 'order');
%
%   See also OPT_MODEL.

%   MP-Opt-Model
%   Copyright (c) 2008-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

val = obj;
for k = 1:length(varargin)
    if ischar(varargin{k})
        val = val.(varargin{k});
    else
        val = val(varargin{k});
    end
end
