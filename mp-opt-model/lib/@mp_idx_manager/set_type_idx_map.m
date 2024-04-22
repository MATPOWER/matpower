function s = set_type_idx_map(obj, set_type, idxs, group_by_name)
% set_type_idx_map - Maps index for set type back to named set & index within set
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

NS = obj.(set_type).NS;

%% special case : everything and grouped by name
if group_by_name && isempty(idxs)
    %% pre-allocate return struct
    c = cell(NS, 1);    %% pre-allocate return cell array
    s = struct('name', c, 'idx', c, 'i', c);
    for k = 1:NS;
        name = obj.(set_type).order(k).name;    %% name of individual set
        idx = obj.(set_type).order(k).idx;      %% idx of individual set
        s(k).name = name;
        if isempty(idx)
            i1 = obj.(set_type).idx.i1.(name);
            iN = obj.(set_type).idx.iN.(name);
            N = obj.(set_type).idx.N.(name);
            s(k).i = [1:N]';
            s(k).j = [i1:iN]';
        else
            s(k).idx = idx;
            ss = substruct('.', name, '()', idx);
            i1 = subsref(obj.(set_type).idx.i1, ss);
            iN = subsref(obj.(set_type).idx.iN, ss);
            N = subsref(obj.(set_type).idx.N, ss);
            s(k).i = [1:N]';
            s(k).j = [i1:iN]';
        end
    end
else    %% general case
    if isempty(idxs)
        idxs = [1:obj.(set_type).N]';   %% all indices
    end

    %% check for invalid idxs
    if any(idxs > obj.(set_type).N)
        error('mp_idx_manager.set_type_idx_map: IDXS must not exceed maximum %s index (%d)', set_type, obj.(set_type).N);
    end
    if any(idxs < 1)
        error('mp_idx_manager.set_type_idx_map: IDXS must be positive');
    end

    %% pre-allocate return struct
    c = cell(size(idxs));       %% pre-allocate return cell array
    s = struct('name', c, 'idx', c, 'i', c);

    %% sort idxs so we can go through loop only once
    [sorted_idxs, jj] = sort(idxs(:), 'descend');
    k = NS;                         %% index into set type (decrementing)
    name = obj.(set_type).order(k).name;    %% name of individual set
    idx = obj.(set_type).order(k).idx;      %% idx of individual set
    for i = 1:length(sorted_idxs)   %% index into sorted_idxs
        ii = sorted_idxs(i);        %% idx of interest
        j = jj(i);                  %% index into s and original idxs
        while k > 0
            if isempty(idx)
                i1 = obj.(set_type).idx.i1.(name);
                if ii >= i1
                    s(j).name = name;
                    s(j).i = ii - i1 + 1;
                    break;
                else
                    k = k - 1;
                    name = obj.(set_type).order(k).name;
                    idx = obj.(set_type).order(k).idx;
                end
            else
                ss = substruct('.', name, '()', idx);
                i1 = subsref(obj.(set_type).idx.i1, ss);
                if ii >= i1
                    s(j).name = name;
                    s(j).i = ii - i1 + 1;
                    s(j).idx = idx;
                    break;
                else
                    k = k - 1;
                    name = obj.(set_type).order(k).name;
                    idx = obj.(set_type).order(k).idx;
                end
            end
        end
    end

    %% consolidate indices into vectors for each unique
    %% name/idx pair, if requested
    if group_by_name
        %% extract fields
        [name, idx, i] = deal(cell(size(idxs)));
        [name{:}] = deal(s.name);
        [idx{:}]  = deal(s.idx);
        [i{:}]    = deal(s.i);      i = cell2mat(i);

        %% find unique name/idx
        name_idx = cellfun(@join_name_idx, name, idx, ...
            'UniformOutput', 0);
        [c, ia, ic] = unique(name_idx);

        %% recreate struct, grouped by name/idx
        c0 = cell(size(c));
        s = struct('name', name(ia), 'idx', idx(ia), 'i', c0, 'j', c0);
        for k = 1:length(ia)
            s(k).i  = i(ic == k);
            s(k).j = idxs(ic == k);
        end
    end
end

function name_idx = join_name_idx(name, idx)
if isempty(idx)
    name_idx = name;
else
    name_idx = [name sprintf('_%d', idx{:})];
end
