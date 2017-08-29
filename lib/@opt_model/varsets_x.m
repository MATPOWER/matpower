function xx = varsets_x(om, x, vs)
%VARSETS_X  Returns a cell array of sub-vectors of X specified by VARSETS
%   X = OM.VARSETS_X(X, VARSETS)
%
%   Returns a cell array of sub-vectors of X specified by VARSETS.
%
%   See also VARSET_LEN

%   MATPOWER
%   Copyright (c) 2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.


if isempty(vs)          %% all rows of x
    xx = x;
else                    %% selected rows of x
    s = struct('type', {'.', '()'}, 'subs', {'', 1});
    xx = cell(size(vs));
    for v = 1:length(vs)
        vname = vs(v).name;
        vidx = vs(v).idx;
        if isempty(vidx)
            j1 = om.var.idx.i1.(vname);     %% starting row in full x
            jN = om.var.idx.iN.(vname);     %% ending row in full x
        else
            % (calls to substruct() are relatively expensive ...
            % s = substruct('.', vname, '()', vidx);
            % ... so replace it with these more efficient lines)
            s(1).subs = vname;
            s(2).subs = vidx;
            j1 = subsref(om.var.idx.i1, s); %% starting row in full x
            jN = subsref(om.var.idx.iN, s); %% ending row in full x
        end
        xx{v} = x(j1:jN);
    end
end
