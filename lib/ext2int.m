function [i2e, bus, gen, branch, areas] = ext2int(bus, gen, branch, areas)
%EXT2INT   Converts external to internal indexing.
%
%   This function has two forms, (1) the old form that operates on
%   and returns individual matrices and (2) the new form that operates
%   on and returns an entire MATPOWER case struct.
%
%   1.  [I2E, BUS, GEN, BRANCH, AREAS] = EXT2INT(BUS, GEN, BRANCH, AREAS)
%       [I2E, BUS, GEN, BRANCH] = EXT2INT(BUS, GEN, BRANCH)
%
%   If the first argument is a matrix, it simply converts from (possibly
%   non-consecutive) external bus numbers to consecutive internal bus
%   numbers which start at 1. Changes are made to BUS, GEN and BRANCH,
%   which are returned along with a vector of indices I2E that can be
%   passed to INT2EXT to perform the reverse conversion, where
%   EXTERNAL_BUS_NUMBER = I2E(INTERNAL_BUS_NUMBER).
%   AREAS is completely ignored and is only included here for backward
%   compatibility of the API.
%
%   Examples:
%       [i2e, bus, gen, branch, areas] = ext2int(bus, gen, branch, areas);
%       [i2e, bus, gen, branch] = ext2int(bus, gen, branch);
%
%   2.  MPC = EXT2INT(MPC)
%       MPC = EXT2INT(MPC, MPOPT)
%
%   If the input is a single MATPOWER case struct, followed optionally
%   by a MATOWER options struct, then all isolated buses, off-line
%   generators and branches are removed along with any generators or
%   branches connected to isolated buses. Then the buses are renumbered
%   consecutively, beginning at 1, and the generators are sorted by
%   increasing bus number. Any 'ext2int' callback routines registered in
%   the case are also invoked automatically. All of the related indexing
%   information and the original data matrices are stored in an 'order'
%   field in the struct to be used by INT2EXT to perform the reverse
%   conversions. If the case is already using internal numbering it is
%   returned unchanged.
%
%   Examples:
%       mpc = ext2int(mpc);
%       mpc = ext2int(mpc, mpopt);
%
%   The 'order' field of MPC used to store the indexing information
%   needed for subsequent internal to external conversion is structured
%   as:
%
%       order
%           state       'i' | 'e'
%           ext | int
%               bus
%               branch
%               gen
%               gencost
%               A
%               N
%           bus
%               e2i
%               i2e
%               status
%                   on
%                   off
%           gen
%               e2i
%               i2e
%               status
%                   on
%                   off
%           branch
%               status
%                   on
%                   off
%
%   See also INT2EXT, E2I_FIELD, E2I_DATA.

%   MATPOWER
%   Copyright (c) 1996-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if isstruct(bus)
    mpc = bus;
    if nargin < 3
        first = ~isfield(mpc, 'order');
        if first || mpc.order.state == 'e'
            %% define names for columns to data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE] = idx_bus;
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS] = idx_gen;
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS] = idx_brch;

            %% initialize order
            if first
                status = struct('on',   [], ...
                                'off',  []  );
                tmp = struct( ...
                        'e2i',      [], ...
                        'i2e',      [], ...
                        'status',   status ...
                    );
                o = struct( ...
                        'ext',      struct( ...
                                'bus',      [], ...
                                'branch',   [], ...
                                'gen',      [] ...
                            ), ...
                        'bus',      tmp, ...
                        'gen',      tmp, ...
                        'branch',   struct('status', status) ...
                    );
            else
                o = mpc.order;
            end

            %% sizes
            nb = size(mpc.bus, 1);
            ng = size(mpc.gen, 1);
            nb0 = nb;
            ng0 = ng;
            if isfield(mpc, 'A') && size(mpc.A, 2) < 2*nb + 2*ng
                dc = 1;
            elseif isfield(mpc, 'N') && size(mpc.N, 2) < 2*nb + 2*ng
                dc = 1;
            else
                dc = 0;
            end

            %% save data matrices with external ordering
            o.ext.bus    = mpc.bus;
            o.ext.branch = mpc.branch;
            o.ext.gen    = mpc.gen;

            %% check that all buses have a valid BUS_TYPE
            bt = mpc.bus(:, BUS_TYPE);
            err = find(~(bt == PQ | bt == PV | bt == REF | bt == NONE));
            if ~isempty(err)
                error('ext2int: bus %d has an invalid BUS_TYPE', err);
            end

            %% determine which buses, branches, gens are connected & in-service
            n2i = sparse(mpc.bus(:, BUS_I), ones(nb, 1), 1:nb, max(mpc.bus(:, BUS_I)), 1);
            bs = (bt ~= NONE);                      %% bus status
            o.bus.status.on     = find(  bs );      %% connected
            o.bus.status.off    = find( ~bs );      %% isolated
            gs = ( mpc.gen(:, GEN_STATUS) > 0 & ... %% gen status
                    bs(n2i(mpc.gen(:, GEN_BUS))) );
            o.gen.status.on     = find(  gs );      %% on and connected
            o.gen.status.off    = find( ~gs );      %% off or isolated
            brs = ( mpc.branch(:, BR_STATUS) & ...  %% branch status
                    bs(n2i(mpc.branch(:, F_BUS))) & ...
                    bs(n2i(mpc.branch(:, T_BUS))) );
            o.branch.status.on  = find(  brs );     %% on and connected
            o.branch.status.off = find( ~brs );

            %% delete stuff that is "out"
            if ~isempty(o.bus.status.off)
                mpc.bus(o.bus.status.off, :) = [];
            end
            if ~isempty(o.branch.status.off)
                mpc.branch(o.branch.status.off, :) = [];
            end
            if ~isempty(o.gen.status.off)
                mpc.gen(o.gen.status.off, :) = [];
            end

            %% update size
            nb = size(mpc.bus, 1);

            %% apply consecutive bus numbering
            o.bus.i2e = mpc.bus(:, BUS_I);
            o.bus.e2i = sparse(o.bus.i2e, 1, 1:nb);
            if nb
                mpc.bus(:, BUS_I)       = o.bus.e2i( mpc.bus(:, BUS_I)      );
                mpc.gen(:, GEN_BUS)     = o.bus.e2i( mpc.gen(:, GEN_BUS)    );
                mpc.branch(:, F_BUS)    = o.bus.e2i( mpc.branch(:, F_BUS)   );
                mpc.branch(:, T_BUS)    = o.bus.e2i( mpc.branch(:, T_BUS)   );
            end

            %% reorder gens in order of increasing bus number
            [tmp, o.gen.i2e] = sort(mpc.gen(:, GEN_BUS));
            [tmp, o.gen.e2i] = sort(o.gen.i2e);
            mpc.gen = mpc.gen(o.gen.i2e, :);

            if isfield(o, 'int')
                o = rmfield(o, 'int');
            end
            o.state = 'i';
            mpc.order = o;

            %% update gencost, bus_name, gentype, genfuel, A and N
            if isfield(mpc, 'gencost')
                ordering = {'gen'};         %% Pg cost only
                if size(mpc.gencost, 1) == 2*ng0
                    ordering{2} = 'gen';    %% include Qg cost
                end
                mpc = e2i_field(mpc, 'gencost', ordering);
            end
            if isfield(mpc, 'bus_name')
                mpc = e2i_field(mpc, 'bus_name', {'bus'});
            end
            if isfield(mpc, 'gentype')
                mpc = e2i_field(mpc, 'gentype', {'gen'});
            end
            if isfield(mpc, 'genfuel')
                mpc = e2i_field(mpc, 'genfuel', {'gen'});
            end
            if isfield(mpc, 'A') || isfield(mpc, 'N')
                if dc
                    ordering = {'bus', 'gen'};
                else
                    ordering = {'bus', 'bus', 'gen', 'gen'};
                end
            end
            if isfield(mpc, 'A')
                mpc = e2i_field(mpc, 'A', ordering, 2);
            end
            if isfield(mpc, 'N')
                mpc = e2i_field(mpc, 'N', ordering, 2);
            end

            %% execute userfcn callbacks for 'ext2int' stage
            if isfield(mpc, 'userfcn')
                if nargin < 2
                    mpopt = struct();
                else
                    mpopt = gen;
                end
                mpc = run_userfcn(mpc.userfcn, 'ext2int', mpc, mpopt);
            end
        end

        i2e = mpc;
    else                    %% convert extra data (DEPRECATED)
        ordering = branch;              %% rename argument
        if nargin < 4
            dim = 1;
        else
            dim = areas;                %% rename argument
        end
        if ischar(gen) || iscell(gen)   %% field
            warning('Calls of the form MPC = EXT2INT(MPC, ''FIELD_NAME'', ...) have been deprecated. Please replace EXT2INT with E2I_FIELD.');
            i2e = e2i_field(mpc, gen, branch, dim);
        else                            %% value
            warning('Calls of the form VAL = EXT2INT(MPC, VAL, ...) have been deprecated. Please replace EXT2INT with E2I_DATA.');
            i2e = e2i_data(mpc, gen, branch, dim);
        end
    end
else            %% old form
    %% define names for columns to data matrices
    [PQ, PV, REF, NONE, BUS_I] = idx_bus;
    [GEN_BUS] = idx_gen;
    [F_BUS, T_BUS] = idx_brch;

    %% create map of external bus numbers to bus indices
    i2e = bus(:, BUS_I);
    e2i = sparse(max(i2e), 1);
    e2i(i2e) = (1:size(bus, 1))';

    %% renumber buses consecutively
    bus(:, BUS_I)               = e2i( bus(:, BUS_I)            );
    gen(:, GEN_BUS)             = e2i( gen(:, GEN_BUS)          );
    branch(:, F_BUS)            = e2i( branch(:, F_BUS)         );
    branch(:, T_BUS)            = e2i( branch(:, T_BUS)         );
end
