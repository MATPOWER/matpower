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
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009-2010 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

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
            s.subs{k} = ':';        %% indexes of all elements in this dimension
        else    %% sA(k) > sB(k)
            s.subs{k} = (1:sB(k));  %% limit indexes to smaller size of B
        end
    end
end

%% do the assignment
A = subsasgn(A, s, B);
