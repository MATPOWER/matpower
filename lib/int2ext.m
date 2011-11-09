function [bus, gen, branch, areas] = int2ext(i2e, bus, gen, branch, areas)
%INT2EXT   Converts internal to external bus numbering.
%
%   This function has two forms, (1) the old form that operates on
%   and returns individual matrices and (2) the new form that operates
%   on and returns an entire MATPOWER case struct.
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
%   See also EXT2INT, I2E_FIELD, I2E_DATA.

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
            error('int2ext: mpc does not have the ''order'' field required for conversion back to external numbering.');
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
            warning('Calls of the form MPC = INT2EXT(MPC, ''FIELD_NAME'', ...) have been deprecated. Please replace INT2EXT with I2E_FIELD.');
            if nargin > 3
                dim = branch;
            else
                dim = 1;
            end
            bus = i2e_field(mpc, bus, gen, dim);
        else                            %% value
            warning('Calls of the form VAL = INT2EXT(MPC, VAL, ...) have been deprecated. Please replace INT2EXT with I2E_DATA.');
            if nargin > 4
                dim = areas;
            else
                dim = 1;
            end
            bus = i2e_data(mpc, bus, gen, branch, dim);
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
