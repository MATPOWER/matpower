function nv = varsets_len(om, vs)
%VARSETS_LEN  Returns the total number of variables in VARSETS
%   NV = OM.VARSETS_LEN(VARSETS)
%
%   Returns the total number of elements in the optimization sub-vector
%   specified by VARSETS.
%
%   See also VARSETS_CELL2STRUCT

%   MATPOWER
%   Copyright (c) 2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

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
        % (calls to substruct() are relatively expensive ...
        % sn = substruct('.', vs(v).name, '()', vs(v).idx);
        % ... so replace it with these more efficient lines)
        sn(1).subs = vs(v).name;
        sn(2).subs = vs(v).idx;
        N = subsref(om.var.idx.N, sn);
        nv = nv + sum(N(:));
    end
end
