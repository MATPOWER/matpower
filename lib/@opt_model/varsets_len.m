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

if isempty(vs)
    nv = om.var.N;
else
    nv = 0;
    s = struct('type', {'.', '()'}, 'subs', {'', 1});
    for v = 1:length(vs)
        % (calls to substruct() are relatively expensive ...
        % s = substruct('.', vs(v).name, '()', vs(v).idx);
        % ... so replace it with these more efficient lines)
        s(1).subs = vs(v).name;
        s(2).subs = vs(v).idx;
        nv = nv + subsref(om.var.idx.N, s);
    end
end
