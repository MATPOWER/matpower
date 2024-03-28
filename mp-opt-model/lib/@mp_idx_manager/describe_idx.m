function label = describe_idx(obj, set_type, idxs)
% describe_idx - Identifies element indices for a give set type.
% ::
%
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
% See also set_type_idx_map, opt_model.

%   MP-Opt-Model
%   Copyright (c) 2012-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

label = cell(size(idxs));       %% pre-allocate return cell array

if ~isempty(idxs)
    s = obj.set_type_idx_map(set_type, idxs);
end
for i = 1:length(idxs(:))
    idx = s(i).idx;
    if isempty(idx)
        label{i} = sprintf('%s(%d)', s(i).name, s(i).i);
    else
        if length(idx) <= 1
            idxstr = sprintf('%d', idx{1});
        else
            idxstr = [sprintf('%d', idx{1}) sprintf(',%d', idx{2:end})];
        end
        label{i} = sprintf('%s(%s)(%d)', s(i).name, idxstr, s(i).i);
    end
end
if isscalar(idxs)               %% return scalar
    label = label{1};
end
