function s = set_type_idx_map(obj, set_type, idxs, group_by_name)
% set_type_idx_map - Maps index for set type back to named set & index within set
%
% .. note::
%    .. deprecated:: 5.0 Please use mp.set_manager.set_type_idx_map instead.
%
% ::
%
%   S = OBJ.SET_TYPE_IDX_MAP(SET_TYPE, IDXS)
%   S = OBJ.SET_TYPE_IDX_MAP(SET_TYPE)
%   S = OBJ.SET_TYPE_IDX_MAP(SET_TYPE, IDXS, GROUP_BY_NAME)
%
%   Returns a struct of same dimensions as IDXS specifying, for each index,
%   the corresponding named set and element within the named set for the
%   specified SET_TYPE. The return struct has the following fields:
%       name : name of corresponding set
%       idx : cell array of indices for the name, if named set is indexed
%       i : index of element within the set
%       j : (only if GROUP_BY_NAME == 1), corresponding index of set type,
%           equal to a particular element of IDXS
%   If IDXS is empty or not provided it defaults to [1:NS]', where NS is the
%   full dimension of the set corresponding to the all elements for the
%   specified set type.
%
%   If GROUP_BY_NAME is true, then the results are consolidated, with a
%   single entry in S for each unique name/idx pair, where the i and j fields
%   are vectors. In this case S is 1 dimensional.
%
%   Examples:
%       s = obj.set_type_idx_map('var', 87));
%       s = obj.set_type_idx_map('lin', [38; 49; 93]));
%       s = obj.set_type_idx_map('var'));
%       s = obj.set_type_idx_map('var', [], 1));
%
%   See also describe_idx, opt_model.

%   MP-Opt-Model
%   Copyright (c) 2012-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% default args
if nargin < 4
    group_by_name = 0;
    if nargin < 3
        idxs = [];
    end
end

s = obj.(set_type).set_type_idx_map(idxs, group_by_name);
