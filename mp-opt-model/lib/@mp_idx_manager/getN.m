function N = getN(obj, set_type, varargin)
%
% .. note::
%    .. deprecated:: 5.0 Please use mp.set_manager.get_N instead.
%
% getN - Returns the number of elements of a given set type.
% ::
%
%   N = OBJ.GETN(SET_TYPE)
%   N = OBJ.GETN(SET_TYPE, NAME)
%   N = OBJ.GETN(SET_TYPE, NAME, IDX_LIST)
%
%   Returns either the total number of elements of a given set type
%   (e.g. variable, constraint, etc.) or the number corresponding to a
%   specified named block, or indexed named set.
%
%   Examples:
%       N = obj.getN('var')      : total number of variables
%       N = obj.getN('lin')      : total number of linear constraints
%       N = obj.getN('var', name)    : # of variables in named set
%       N = obj.getN('lin', name)    : # of linear constraints in named set
%       N = obj.getN('var', name, idx_list) : # of variables in indexed named set
%
% See also opt_model.

%   MP-Opt-Model
%   Copyright (c) 2008-2025, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

N = obj.(set_type).get_N(varargin{:});
