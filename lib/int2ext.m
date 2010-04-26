function [bus, gen, branch, areas] = int2ext(i2e, bus, gen, branch, areas)
%INT2EXT   Converts internal to external bus numbering.
%
%   This function performs several different tasks, depending on the
%   arguments passed.
%
%   1.  [BUS, GEN, BRANCH, AREAS] = INT2EXT(I2E, BUS, GEN, BRANCH, AREAS)
%       [BUS, GEN, BRANCH] = INT2EXT(I2E, BUS, GEN, BRANCH)
%
%   Converts from the consecutive internal bus numbers back to the originals
%   using the mapping provided by the I2E vector returned from EXT2INT,
%   where EXTERNAL_BUS_NUMBER = I2E(INTERNAL_BUS_NUMBER).
%
%   Examples:
%       [bus, gen, branch, areas] = int2ext(i2e, bus, gen, branch, areas);
%       [bus, gen, branch] = int2ext(i2e, bus, gen, branch);
%
%   2.  MPC = INT2EXT(MPC)
%
%   If the input is a single MATPOWER case struct, then it restores all
%   buses, generators and branches that were removed because of being
%   isolated or off-line, and reverts to the original generator ordering
%   and original bus numbering. This requires that the 'order' field
%   created by EXT2INT be in place.
%
%   Example:
%       mpc = int2ext(mpc);
%
%   3.  VAL = INT2EXT(MPC, VAL, OLDVAL, ORDERING)
%       VAL = INT2EXT(MPC, VAL, OLDVAL, ORDERING, DIM)
%       MPC = INT2EXT(MPC, FIELD, ORDERING)
%       MPC = INT2EXT(MPC, FIELD, ORDERING, DIM)
%
%   For a case struct using internal indexing, this function can be
%   used to convert other data structures as well by passing in 2 to 4
%   extra parameters in addition to the case struct. If the values passed
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
%   If the 2nd argument is a string or cell array of strings, it
%   specifies a field in the case struct whose value should be
%   converted as described above. In this case, the corresponding
%   OLDVAL is taken from where it was stored by EXT2INT in
%   MPC.ORDER.EXT and the updated case struct is returned.
%   If FIELD is a cell array of strings, they specify nested fields.
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
%       A_ext = int2ext(mpc, A_int, A_orig, {'bus','bus','gen','gen'}, 2);
%
%       Converts an A matrix for user-supplied OPF constraints from
%       internal to external ordering, where the columns of the A
%       matrix correspond to bus voltage angles, then voltage
%       magnitudes, then generator real power injections and finally
%       generator reactive power injections.
%
%       gencost_ext = int2ext(mpc, gencost_int, gencost_orig, {'gen','gen'}, 1);
%
%       Converts a GENCOST matrix that has both real and reactive power
%       costs (in rows 1--ng and ng+1--2*ng, respectively).
%
%       mpc = int2ext(mpc, {'reserves', 'cost'}, 'gen');
%
%       Reorders rows of mpc.reserves.cost to match external generator
%       ordering.
%
%       mpc = int2ext(mpc, {'reserves', 'zones'}, 'gen', 2);
%
%       Reorders columns of mpc.reserves.zones to match external
%       generator ordering.
%
%   See also EXT2INT.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2010 by Power System Engineering Research Center (PSERC)
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

if isstruct(i2e)
    mpc = i2e;
    if nargin == 1
        if ~isfield(mpc, 'order')
            error('int2ext: mpc does not have the ''order'' field require for conversion back to external numbering.');
        end
        o = mpc.order;

        if o.state == 'i'
            %% define names for columns to data matrices
            [PQ, PV, REF, NONE, BUS_I] = idx_bus;
            GEN_BUS = idx_gen;
            [F_BUS, T_BUS] = idx_brch;
            [AREA_I, PRICE_REF_BUS] = idx_area;

            %% execute userfcn callbacks for 'int2ext' stage
            if isfield(mpc, 'userfcn')
                mpc = run_userfcn(mpc.userfcn, 'int2ext', mpc);
            end

            %% save data matrices with internal ordering & restore originals
            o.int.bus    = mpc.bus;
            o.int.branch = mpc.branch;
            o.int.gen    = mpc.gen;
            mpc.bus     = o.ext.bus;
            mpc.branch  = o.ext.branch;
            mpc.gen     = o.ext.gen;
            if isfield(mpc, 'gencost')
                o.int.gencost = mpc.gencost;
                mpc.gencost = o.ext.gencost;
            end
            if isfield(mpc, 'areas')
                o.int.areas = mpc.areas;
                mpc.areas = o.ext.areas;
            end
            if isfield(mpc, 'A')
                o.int.A = mpc.A;
                mpc.A = o.ext.A;
            end
            if isfield(mpc, 'N')
                o.int.N = mpc.N;
                mpc.N = o.ext.N;
            end

            %% update data (in bus, branch and gen only)
            mpc.bus(o.bus.status.on, :)       = o.int.bus;
            mpc.branch(o.branch.status.on, :) = o.int.branch;
            mpc.gen(o.gen.status.on, :)       = o.int.gen(o.gen.i2e, :);
            if isfield(mpc, 'areas')
                mpc.areas(o.areas.status.on, :) = o.int.areas;
            end

            %% revert to original bus numbers
            mpc.bus(o.bus.status.on, BUS_I) = ...
                    o.bus.i2e( mpc.bus(o.bus.status.on, BUS_I) );
            mpc.branch(o.branch.status.on, F_BUS) = ...
                    o.bus.i2e( mpc.branch(o.branch.status.on, F_BUS) );
            mpc.branch(o.branch.status.on, T_BUS) = ...
                    o.bus.i2e( mpc.branch(o.branch.status.on, T_BUS) );
            mpc.gen(o.gen.status.on, GEN_BUS) = ...
                    o.bus.i2e( mpc.gen(o.gen.status.on, GEN_BUS) );
            if isfield(mpc, 'areas')
                mpc.areas(o.areas.status.on, PRICE_REF_BUS) = ...
                        o.bus.i2e( mpc.areas(o.areas.status.on, PRICE_REF_BUS) );
            end

            if isfield(o, 'ext')
                o = rmfield(o, 'ext');
            end
            o.state = 'e';
            mpc.order = o;
        else
            error('int2ext: mpc claims it is already using external numbering.');
        end

        bus = mpc;
    else                    %% convert extra data
        if ischar(bus) || iscell(bus)   %% field
            field = bus;
            ordering = gen;
            if nargin > 3
                dim = branch;
            else
                dim = 1;
            end
            if ischar(field)
                mpc.order.int.(field) = mpc.(field);
                mpc.(field) = int2ext(mpc, mpc.(field), ...
                                mpc.order.ext.(field), ordering, dim);
            else
                for k = 1:length(field)
                    s(k).type = '.';
                    s(k).subs = field{k};
                end
                if ~isfield(mpc.order, 'int')
                    mpc.order.int = [];
                end
                mpc.order.int = subsasgn(mpc.order.int, s, subsref(mpc, s));
                mpc = subsasgn(mpc, s, int2ext(mpc, subsref(mpc, s), ...
                    subsref(mpc.order.ext, s), ordering, dim));
            end
            bus = mpc;
        else                            %% value
            val = bus;
            oldval = gen;
            ordering = branch;
            o = mpc.order;
            if nargin > 4
                dim = areas;
            else
                dim = 1;
            end
            if ischar(ordering)         %% single set
                if strcmp(ordering, 'gen')
                    v = get_reorder(val, o.(ordering).i2e, dim);
                else
                    v = val;
                end
                bus = set_reorder(oldval, v, o.(ordering).status.on, dim);
            else                            %% multiple sets
                be = 0;  %% base, external indexing
                bi = 0;  %% base, internal indexing
                for k = 1:length(ordering)
                    ne = size(o.ext.(ordering{k}), 1);
                    ni = size(mpc.(ordering{k}), 1);
                    v = get_reorder(val, bi+(1:ni), dim);
                    oldv = get_reorder(oldval, be+(1:ne), dim);
                    new_v{k} = int2ext(mpc, v, oldv, ordering{k}, dim);
                    be = be + ne;
                    bi = bi + ni;
                end
                ni = size(val, dim);
                if ni > bi              %% the rest
                    v = get_reorder(val, bi+1:ni, dim);
                    new_v{length(new_v)+1} = v;
                end
                bus = cat(dim, new_v{:});
            end
        end
    end
else            %% old form
    %% define names for columns to data matrices
    [PQ, PV, REF, NONE, BUS_I] = idx_bus;
    [GEN_BUS] = idx_gen;
    [F_BUS, T_BUS] = idx_brch;
    [AREA_I, PRICE_REF_BUS] = idx_area;

    bus(:, BUS_I)               = i2e( bus(:, BUS_I)            );
    gen(:, GEN_BUS)             = i2e( gen(:, GEN_BUS)          );
    branch(:, F_BUS)            = i2e( branch(:, F_BUS)         );
    branch(:, T_BUS)            = i2e( branch(:, T_BUS)         );
    if nargin > 4 && nargout > 3 && ~isempty(areas)
        areas(:, PRICE_REF_BUS) = i2e( areas(:, PRICE_REF_BUS)  );
    end
end
