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
%   The 2nd argument is a string or cell array of strings, specifying
%   a field in the case struct whose value should be converted by
%   a corresponding call to I2E_DATA. The field can contain either a
%   numeric or a cell array. The corresponding OLDVAL is taken from
%   where it was stored by EXT2INT in MPC.ORDER.EXT and the updated
%   case struct is returned. If FIELD is a cell array of strings,
%   they specify nested fields.
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
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

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
