function mpc = e2i_field(mpc, field, ordering, dim)
%E2I_FIELD   Converts fields of MPC from external to internal indexing.
%
%   This function performs several different tasks, depending on the
%   arguments passed.
%
%   MPC = E2I_FIELD(MPC, FIELD, ORDERING)
%   MPC = E2I_FIELD(MPC, FIELD, ORDERING, DIM)
%
%   When given a case struct that has already been converted to
%   internal indexing, this function can be used to convert other data
%   structures as well by passing in 2 or 3 extra parameters in
%   addition to the case struct.
%
%   The 2nd argument is a string or cell array of strings, specifying
%   a field in the case struct whose value should be converted by
%   a corresponding call to E2I_DATA. The field can contain either a
%   numeric or a cell array. The converted value is stored back in the
%   specified field, the original value is saved for later use and the
%   updated case struct is returned. If FIELD is a cell array of strings,
%   they specify nested fields.
%
%   The 3rd and optional 4th arguments are simply passed along to
%   the call to E2I_DATA.
%
%   Examples:
%       mpc = e2i_field(mpc, {'reserves', 'cost'}, 'gen');
%
%       Reorders rows of mpc.reserves.cost to match internal generator
%       ordering.
%
%       mpc = e2i_field(mpc, {'reserves', 'zones'}, 'gen', 2);
%
%       Reorders columns of mpc.reserves.zones to match internal
%       generator ordering.
%
%   See also I2E_FIELD, E2I_DATA, EXT2INT.

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
    mpc.order.ext.(field) = mpc.(field);
    mpc.(field) = e2i_data(mpc, mpc.(field), ordering, dim);
else    %% iscell(field)
    for k = 1:length(field)
        s(k).type = '.';
        s(k).subs = field{k};
    end
    mpc.order.ext = subsasgn(mpc.order.ext, s, subsref(mpc, s));
    mpc = subsasgn(mpc, s, ...
        e2i_data(mpc, subsref(mpc, s), ordering, dim) );
end
