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
%   AREAS is completely ignored and is only included here for backward
%   compatibility of the API.
%
%   Examples:
%       [bus, gen, branch, areas] = int2ext(i2e, bus, gen, branch, areas);
%       [bus, gen, branch] = int2ext(i2e, bus, gen, branch);
%
%   2.  MPC = INT2EXT(MPC)
%       MPC = INT2EXT(MPC, MPOPT)
%
%   If the input is a single MATPOWER case struct, followed optionally
%   by a MATOWER options struct, then it restores all buses, generators
%   and branches that were removed because of being isolated or off-line,
%   and reverts to the original generator ordering and original bus
%   numbering. This requires that the 'order' field created by EXT2INT be
%   in place.
%
%   Examples:
%       mpc = int2ext(mpc);
%       mpc = int2ext(mpc, mpopt);
%
%   See also EXT2INT, I2E_FIELD, I2E_DATA.

%   MATPOWER
%   Copyright (c) 1996-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if isstruct(i2e)
    mpc = i2e;
    if nargin < 3
        if ~isfield(mpc, 'order')
            error('int2ext: mpc does not have the ''order'' field required for conversion back to external numbering.');
        end
        o = mpc.order;

        if o.state == 'i'
            %% define names for columns to data matrices
            [PQ, PV, REF, NONE, BUS_I] = idx_bus;
            GEN_BUS = idx_gen;
            [F_BUS, T_BUS] = idx_brch;

            %% execute userfcn callbacks for 'int2ext' stage
            if isfield(mpc, 'userfcn')
                if nargin < 2
                    mpopt = struct();
                else
                    mpopt = gen;
                end
                mpc = run_userfcn(mpc.userfcn, 'int2ext', mpc, mpopt);
            end
            
            %% convert back "extra" fields
            if isfield(mpc, 'gencost')
                ordering = {'gen'};         %% Pg cost only
                if size(mpc.gencost, 1) == 2*size(mpc.gen, 1) && ...
                        size(mpc.gencost, 1) ~= 0
                    ordering{2} = 'gen';    %% include Qg cost
                end
                mpc = i2e_field(mpc, 'gencost', ordering);
            end
            if isfield(mpc, 'bus_name')
                mpc = i2e_field(mpc, 'bus_name', {'bus'});
            end
            if isfield(mpc, 'gentype')
                mpc = i2e_field(mpc, 'gentype', {'gen'});
            end
            if isfield(mpc, 'genfuel')
                mpc = i2e_field(mpc, 'genfuel', {'gen'});
            end
            %% assume A and N are "read-only"
            %% (otherwise need to convert back, using i2e_field() which
            %% requires knowing if they are sized for AC or DC)
            if isfield(mpc, 'A')
                o.int.A = mpc.A;
                mpc.A = o.ext.A;
            end
            if isfield(mpc, 'N')
                o.int.N = mpc.N;
                mpc.N = o.ext.N;
            end

            %% save data matrices with internal ordering & restore originals
            o.int.bus    = mpc.bus;
            o.int.branch = mpc.branch;
            o.int.gen    = mpc.gen;
            mpc.bus     = o.ext.bus;
            mpc.branch  = o.ext.branch;
            mpc.gen     = o.ext.gen;

            %% zero pad data matrices on right if necessary
            nci = size(o.int.bus, 2);
            [nr, nc] = size(mpc.bus);
            if nc < nci
                mpc.bus = [mpc.bus zeros(nr, nci-nc)];
            end
            nci = size(o.int.branch, 2);
            [nr, nc] = size(mpc.branch);
            if nc < nci
                mpc.branch = [mpc.branch zeros(nr, nci-nc)];
            end
            nci = size(o.int.gen, 2);
            [nr, nc] = size(mpc.gen);
            if nc < nci
                mpc.gen = [mpc.gen zeros(nr, nci-nc)];
            end

            %% update data (in bus, branch, and gen only)
            mpc.bus(o.bus.status.on, :)       = o.int.bus;
            mpc.branch(o.branch.status.on, :) = o.int.branch;
            mpc.gen(o.gen.status.on, :)       = o.int.gen(o.gen.e2i, :);

            %% revert to original bus numbers
            mpc.bus(o.bus.status.on, BUS_I) = ...
                    o.bus.i2e( mpc.bus(o.bus.status.on, BUS_I) );
            mpc.branch(o.branch.status.on, F_BUS) = ...
                    o.bus.i2e( mpc.branch(o.branch.status.on, F_BUS) );
            mpc.branch(o.branch.status.on, T_BUS) = ...
                    o.bus.i2e( mpc.branch(o.branch.status.on, T_BUS) );
            mpc.gen(o.gen.status.on, GEN_BUS) = ...
                    o.bus.i2e( mpc.gen(o.gen.status.on, GEN_BUS) );

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

    bus(:, BUS_I)               = i2e( bus(:, BUS_I)            );
    gen(:, GEN_BUS)             = i2e( gen(:, GEN_BUS)          );
    branch(:, F_BUS)            = i2e( branch(:, F_BUS)         );
    branch(:, T_BUS)            = i2e( branch(:, T_BUS)         );
end
