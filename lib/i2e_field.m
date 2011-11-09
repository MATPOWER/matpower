function mpc = i2e_field(mpc, field, ordering, dim)
%I2E_FIELD   Converts fields of MPC from internal to external bus numbering.
%
%   MPC = I2E_FIELD(MPC, FIELD, ORDERING)
%   MPC = I2E_FIELD(MPC, FIELD, ORDERING, DIM)
%
%   For a case struct using internal indexing, this function can be
%   used to convert other data structures as well by passing in 2 or 3
%   extra parameters in addition to the case struct.
%
%   If the 2nd argument is a string or cell array of strings, it
%   specifies a field in the case struct whose value should be
%   converted by I2E_DATA. In this case, the corresponding
%   OLDVAL is taken from where it was stored by EXT2INT in
%   MPC.ORDER.EXT and the updated case struct is returned.
%   If FIELD is a cell array of strings, they specify nested fields.
%
%   The 3rd and optional 4th arguments are simply passed along to
%   the call to I2E_DATA.
%
%   Examples:
%       mpc = i2e_field(mpc, {'reserves', 'cost'}, 'gen');
%
%       Reorders rows of mpc.reserves.cost to match external generator
%       ordering.
%
%       mpc = i2e_field(mpc, {'reserves', 'zones'}, 'gen', 2);
%
%       Reorders columns of mpc.reserves.zones to match external
%       generator ordering.
%
%   See also E2I_FIELD, I2E_DATA, INT2EXT.

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

if nargin < 4
    dim = 1;
end
if ischar(field)
    mpc.order.int.(field) = mpc.(field);
    mpc.(field) = i2e_data(mpc, mpc.(field), ...
                    mpc.order.ext.(field), ordering, dim);
else    %% iscell(field)
    for k = 1:length(field)
        s(k).type = '.';
        s(k).subs = field{k};
    end
    if ~isfield(mpc.order, 'int')
        mpc.order.int = [];
    end
    mpc.order.int = subsasgn(mpc.order.int, s, subsref(mpc, s));
    mpc = subsasgn(mpc, s, i2e_data(mpc, subsref(mpc, s), ...
        subsref(mpc.order.ext, s), ordering, dim));
end
