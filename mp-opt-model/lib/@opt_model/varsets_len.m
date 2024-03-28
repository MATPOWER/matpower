function nv = varsets_len(om, vs)
% varsets_len - Returns the total number of variables in VARSETS
% ::
%
%   NV = OM.VARSETS_LEN(VARSETS)
%
%   Returns the total number of elements in the optimization sub-vector
%   specified by VARSETS.
%
% See also varsets_cell2struct.

%   MP-Opt-Model
%   Copyright (c) 2017-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

persistent sn;
if isempty(vs)
    nv = om.var.N;
else
    nv = 0;

    %% calls to substruct() are relatively expensive, so we pre-build the
    %% struct for addressing numeric array fields, updating only
    %% the subscripts before use
    if isempty(sn)
        sn = struct('type', {'.', '()'}, 'subs', {'', 1});
    end

    for v = 1:length(vs)
        idx = vs(v).idx;
        if isempty(idx)
            N = om.var.idx.N.(vs(v).name);
        else
            % (calls to substruct() are relatively expensive ...
            % sn = substruct('.', vs(v).name, '()', vs(v).idx);
            % ... so replace it with these more efficient lines)
            sn(1).subs = vs(v).name;
            sn(2).subs = idx;
            N = subsref(om.var.idx.N, sn);
        end
        nv = nv + sum(N(:));
    end
end
