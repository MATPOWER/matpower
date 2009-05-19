function A = set_reorder(A, B, idx, dim)
%SET_REORDER Assigns B to A with one of the dimensions of A indexed.
%
%   A = set_reorder(A, B, idx, dim)
%
%   Returns A after doing A(:, ..., :, idx, :, ..., :) = B
%   where dim determines in which dimension to place the idx.

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
A = subsasgn(A, s, B);

return;
