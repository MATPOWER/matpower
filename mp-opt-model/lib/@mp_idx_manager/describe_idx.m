function label = describe_idx(obj, set_type, idxs)
%DESCRIBE_IDX  Identifies element indices for a give set type.
%   LABEL = OBJ.DESCRIBE_IDX(SET_TYPE, IDXS)
%
%   Returns strings describing (name and index) the element of the
%   specified set type (e.g. variable or constraint) that corresponds to
%   the indices in IDXS. The return value is a string if IDXS is a scalar,
%   otherwise it is a cell array of strings of the same dimension as IDXS.
%
%   Examples:
%       label = obj.describe_idx('var', 87));
%       labels = obj.describe_idx('lin', [38; 49; 93]));
%
%   See also OPT_MODEL.

%   MP-Opt-Model
%   Copyright (c) 2012-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

label = cell(size(idxs));       %% pre-allocate return cell array
for i = 1:length(idxs(:))
    ii = idxs(i);
    if ii > obj.(set_type).N
        error('@mp_idx_manager/describe_idx: index exceeds maximum %s index (%d)', set_type, obj.(set_type).N);
    end
    if ii < 1
        error('@mp_idx_manager/describe_idx: index must be positive');
    end
    for k = obj.(set_type).NS:-1:1
        name = obj.(set_type).order(k).name;
        idx = obj.(set_type).order(k).idx;
        if isempty(idx)
            if ii >= obj.(set_type).idx.i1.(name)
                label{i} = sprintf('%s(%d)', name, ii - obj.(set_type).idx.i1.(name) + 1);
                break;
            end
        else
            s = substruct('.', name, '()', idx);
            if ii >= subsref(obj.(set_type).idx.i1, s)
                idxstr = sprintf('%d', idx{1});
                for j = 2:length(idx)
                    idxstr = sprintf('%s,%d', idxstr, idx{j});
                end
                label{i} = sprintf('%s(%s)(%d)', name, idxstr, ...
                            ii - subsref(obj.(set_type).idx.i1, s) + 1);
                break;
            end
        end
    end
end
if isscalar(idxs)               %% return scalar
    label = label{1};
end
