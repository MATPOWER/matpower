function xx = varsets_x(om, x, vs, want_vector)
%VARSETS_X  Returns a cell array of sub-vectors of X specified by VARSETS
%   X = OM.VARSETS_X(X, VARSETS)
%   X = OM.VARSETS_X(X, VARSETS, 'vector')
%
%   Returns a cell array of sub-vectors of X specified by VARSETS, or
%   the full optimization vector X, if VARSETS is empty.
%
%   If a 3rd argument is present (value is ignored) the returned value is
%   a single numeric vector.
%
%   See also VARSET_LEN

%   MATPOWER
%   Copyright (c) 2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.


persistent sn;
if isempty(vs)          %% all rows of x
    xx = x;
else                    %% selected rows of x
    vsN = length(vs);
    xx = cell(1, vsN);

    %% calls to substruct() are relatively expensive, so we pre-build the
    %% struct for addressing numeric array fields, updating only
    %% the subscripts before use
    if isempty(sn)
        sn = struct('type', {'.', '()'}, 'subs', {'', 1});
    end

    ki = 0;
    for v = 1:length(vs)
        vname = vs(v).name;
        vidx = vs(v).idx;
        if isempty(vidx)
            i1 = om.var.idx.i1.(vname);         %% starting row in full x
            iN = om.var.idx.iN.(vname);         %% ending row in full x
        else
            % (calls to substruct() are relatively expensive ...
            % sn = substruct('.', vname, '()', vidx);
            % ... so replace it with these more efficient lines)
            sn(1).subs = vname;
            sn(2).subs = vidx;
            i1 = subsref(om.var.idx.i1, sn);    %% starting row in full x
            iN = subsref(om.var.idx.iN, sn);    %% ending row in full x
        end
        if isscalar(i1)             %% simple named set, or indexed named set
            ki = ki + 1;
            xx{ki} = x(i1:iN);              %% single set of indices for varset
        else                        %% multi-dim named set w/no index specified
            ii1 = permute(i1, ndims(i1):-1:1);
            iiN = permute(iN, ndims(i1):-1:1);
            for j = 1:length(ii1(:))
                ki = ki + 1;
                xx{ki} = x(ii1(j):iiN(j));  %% multiple sets of indices for varset
            end
        end
    end
    if nargin > 3
        xx = vertcat(xx{:});
    end
end
