function label = describe_idx(om, set_type, idxs)
%DESCRIBE_IDX  Identifies variable, constraint and cost row indices.
%   LABEL = OM.DESCRIBE_IDX(SET_TYPE, IDXS)
%
%   Returns strings describing (name and index) the variable, constraint
%   or cost row that corresponds to the indices in IDXS. SET_TYPE must be
%   one of the following: 'var', 'lin', 'nle', 'nli' or 'cost',
%   corresponding to indices for variables, linear constraints, nonlinear
%   equality constraints, nonlinear inequality constraints and cost rows,
%   respectively. The return value is a string if IDXS is a scalar,
%   otherwise it is a cell array of strings of the same dimension as IDXS.
%
%   Examples:
%       label = om.describe_idx('var', 87));
%       labels = om.describe_idx('lin', [38; 49; 93]));
%
%   See also OPT_MODEL.

%   MATPOWER
%   Copyright (c) 2012-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

label = cell(size(idxs));       %% pre-allocate return cell array
for i = 1:length(idxs(:))
    ii = idxs(i);
    if ii > om.(set_type).N
        error('@opt_model/describe_idx: index exceeds maximum %s index (%d)', set_type, om.(set_type).N);
    end
    if ii < 1
        error('@opt_model/describe_idx: index must be positive');
    end
    for k = om.(set_type).NS:-1:1
        name = om.(set_type).order(k).name;
        idx = om.(set_type).order(k).idx;
        if isempty(idx)
            if ii >= om.(set_type).idx.i1.(name)
                label{i} = sprintf('%s(%d)', name, ii - om.(set_type).idx.i1.(name) + 1);
                break;
            end
        else
            s = substruct('.', name, '()', idx);
            if ii >= subsref(om.(set_type).idx.i1, s)
                idxstr = sprintf('%d', idx{1});
                for j = 2:length(idx)
                    idxstr = sprintf('%s,%d', idxstr, idx{j});
                end
                label{i} = sprintf('%s(%s)(%d)', name, idxstr, ...
                            ii - subsref(om.(set_type).idx.i1, s) + 1);
                break;
            end
        end
    end
end
if isscalar(idxs)               %% return scalar
    label = label{1};
end
