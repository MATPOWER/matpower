function A = set_reorder(A, B, idx, dim)
%SET_REORDER Assigns B to A with one of the dimensions of A indexed.
%
%   A = SET_REORDER(A, B, IDX, DIM)
%
%   Returns A after doing A(:, ..., :, IDX, :, ..., :) = B
%   where DIM determines in which dimension to place the IDX.
%
%   See also GET_REORDER.

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
