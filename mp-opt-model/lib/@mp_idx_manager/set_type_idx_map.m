function s = set_type_idx_map(obj, set_type, idxs)
%SET_TYPE_IDX_MAP  Maps index for set type back to named set & index within set
%   S = OBJ.SET_TYPE_IDX_MAP(SET_TYPE, IDXS)
%
%   Returns a struct of same dimensions as IDXS specifying, for each index,
%   the corresponding named set and element within the named set for the
%   specified SET_TYPE. The return struct has the following fields:
%       name : name of corresponding set
%       idx : cell array of indices for the name, if named set is indexed
%       i : index of element within the set
%
%   Examples:
%       s = obj.set_type_idx_map('var', 87));
%       s = obj.set_type_idx_map('lin', [38; 49; 93]));
%
%   See also DESCRIBE_IDX, OPT_MODEL.

%   MP-Opt-Model
%   Copyright (c) 2012-2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% pre-allocate return struct
c = cell(size(idxs));       %% pre-allocate return cell array
s = struct('name', c, 'idx', c, 'i', c);

for i = 1:length(idxs(:))
    ii = idxs(i);
    if ii > obj.(set_type).N
        error('@mp_idx_manager/set_type_idx_map: index exceeds maximum %s index (%d)', set_type, obj.(set_type).N);
    end
    if ii < 1
        error('@mp_idx_manager/set_type_idx_map: index must be positive');
    end
    for k = obj.(set_type).NS:-1:1
        name = obj.(set_type).order(k).name;
        idx = obj.(set_type).order(k).idx;
        if isempty(idx)
            if ii >= obj.(set_type).idx.i1.(name)
                s(i).name = name;
                s(i).i = ii - obj.(set_type).idx.i1.(name) + 1;
                break;
            end
        else
            ss = substruct('.', name, '()', idx);
            if ii >= subsref(obj.(set_type).idx.i1, ss)
                s(i).name = name;
                s(i).i = ii - subsref(obj.(set_type).idx.i1, ss) + 1;
                s(i).idx = idx;
                break;
            end
        end
    end
end
