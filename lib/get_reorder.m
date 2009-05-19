function B = get_reorder(A, idx, dim)
%GET_REORDER    Returns A with one of its dimensions indexed.
%
%   B = get_reorder(A, idx, dim)
%
%   Returns A(:, ..., :, idx, :, ..., :), where dim determines
%   in which dimension to place the idx.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009 by Power System Engineering Research Center (PSERC)
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

return;
