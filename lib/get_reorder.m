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
%   other modules (such as M-files and MEX-files) available in a
%   Matlab (or compatible) environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

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
