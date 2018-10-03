function A = set_reorder(A, B, idx, dim)
%SET_REORDER Assigns B to A with one of the dimensions of A indexed.
%
%   A = SET_REORDER(A, B, IDX, DIM)
%
%   Returns A after doing A(:, ..., :, IDX, :, ..., :) = B
%   where DIM determines in which dimension to place the IDX.
%
%   If any dimension of B is smaller than the corresponding dimension
%   of A, the "extra" elements in A are untouched. If any dimension of
%   B is larger than the corresponding dimension of A, then A is padded
%   with zeros (if numeric) or empty matrices (if cell array) before
%   performing the assignment.
%
%   See also GET_REORDER.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% check dimensions
ndim = ndims(A);
sA = size(A);
sB = size(B);
d = (1:length(sA));
d(dim) = [];        %% indices of all dimensions other than DIM

%% pad A with zeros (numeric) or empty matrices (cell), if necessary
s.subs = cell(1, ndim);
if any(sA(d) < sB(d))
    s.subs = num2cell(max(sA, sB));
    if iscell(A)
        s.type = '{}';
        A = subsasgn(A, s, []);
    else
        s.type = '()';
        A = subsasgn(A, s, 0);
    end
end

%% set up indexing
s.type = '()';
for k = 1:ndim
    if k == dim
        s.subs{k} = idx;
    else
        if sA(k) == sB(k)
            s.subs{k} = ':';        %% indices of all elements in this dimension
        else    %% sA(k) > sB(k)
            s.subs{k} = (1:sB(k));  %% limit indices to smaller size of B
        end
    end
end

%% do the assignment
A = subsasgn(A, s, B);
