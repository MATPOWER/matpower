function newval = i2e_data(mpc, val, oldval, ordering, dim)
%I2E_DATA   Converts data from internal to external bus numbering.
%
%   VAL = I2E_DATA(MPC, VAL, OLDVAL, ORDERING)
%   VAL = I2E_DATA(MPC, VAL, OLDVAL, ORDERING, DIM)
%
%   For a case struct using internal indexing, this function can be
%   used to convert other data structures as well by passing in 3 or 4
%   extra parameters in addition to the case struct. If the value passed
%   in the 2nd argument (VAL) is a column vector, it will be converted
%   according to the ordering specified by the 4th argument (ORDERING,
%   described below). If VAL is an n-dimensional matrix, then the
%   optional 5th argument (DIM, default = 1) can be used to specify
%   which dimension to reorder. The 3rd argument (OLDVAL) is used to
%   initialize the return value before converting VAL to external
%   indexing. In particular, any data corresponding to off-line gens
%   or branches or isolated buses or any connected gens or branches
%   will be taken from OLDVAL, with VAL supplying the rest of the
%   returned data.
%
%   The ORDERING argument is used to indicate whether the data
%   corresponds to bus-, gen- or branch-ordered data. It can be one
%   of the following three strings: 'bus', 'gen' or 'branch'. For
%   data structures with multiple blocks of data, ordered by bus,
%   gen or branch, they can be converted with a single call by
%   specifying ORDERING as a cell array of strings.
%
%   Any extra elements, rows, columns, etc. beyond those indicated
%   in ORDERING, are not disturbed.
%
%   Examples:
%       A_ext = i2e_data(mpc, A_int, A_orig, {'bus','bus','gen','gen'}, 2);
%
%       Converts an A matrix for user-supplied OPF constraints from
%       internal to external ordering, where the columns of the A
%       matrix correspond to bus voltage angles, then voltage
%       magnitudes, then generator real power injections and finally
%       generator reactive power injections.
%
%       gencost_ext = i2e_data(mpc, gencost_int, gencost_orig, {'gen','gen'}, 1);
%
%       Converts a GENCOST matrix that has both real and reactive power
%       costs (in rows 1--ng and ng+1--2*ng, respectively).
%
%   See also E2I_DATA, I2E_FIELD, INT2EXT.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009-2011 by Power System Engineering Research Center (PSERC)
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

if ~isfield(mpc, 'order')
    error('i2e_data: mpc does not have the ''order'' field required for conversion back to external numbering.');
end
o = mpc.order;
if o.state ~= 'i'
    error('i2e_data: mpc does not appear to be in internal order');
end
if nargin < 5
    dim = 1;
end
if ischar(ordering)         %% single set
    if strcmp(ordering, 'gen')
        v = get_reorder(val, o.(ordering).i2e, dim);
    else
        v = val;
    end
    newval = set_reorder(oldval, v, o.(ordering).status.on, dim);
else                        %% multiple sets
    be = 0;  %% base, external indexing
    bi = 0;  %% base, internal indexing
    for k = 1:length(ordering)
        ne = size(o.ext.(ordering{k}), 1);
        ni = size(mpc.(ordering{k}), 1);
        v = get_reorder(val, bi+(1:ni), dim);
        oldv = get_reorder(oldval, be+(1:ne), dim);
        new_v{k} = i2e_data(mpc, v, oldv, ordering{k}, dim);
        be = be + ne;
        bi = bi + ni;
    end
    ni = size(val, dim);
    if ni > bi              %% the rest
        v = get_reorder(val, bi+1:ni, dim);
        new_v{length(new_v)+1} = v;
    end
    newval = cat(dim, new_v{:});
end
