function B = get_reorder(A, idx, dim)
%GET_REORDER    Returns A with one of its dimensions indexed.
%
%   B = GET_REORDER(A, IDX, DIM)
%
%   Returns A(:, ..., :, IDX, :, ..., :), where DIM determines
%   in which dimension to place the IDX.
%
%   See also SET_REORDER.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

ndim = ndims(A);
s.type = '()';
s.subs = cell(1, ndim);
for k = 1:ndim
    if k == dim
        s.subs{k} = idx;
    else
        s.subs{k} = ':';
    end
end
B = subsref(A, s);
