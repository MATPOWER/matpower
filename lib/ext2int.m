function [i2e, bus, gen, branch, areas] = ext2int(bus, gen, branch, areas)
%EXT2INT   Converts external to internal indexing.
%
%   This function performs several different tasks, depending on the
%   arguments passed.
%
%   1.  [I2E, BUS, GEN, BRANCH, AREAS] = EXT2INT(BUS, GEN, BRANCH, AREAS)
%       [I2E, BUS, GEN, BRANCH] = EXT2INT(BUS, GEN, BRANCH)
%
%   If the first argument is a matrix, it simply converts from (possibly
%   non-consecutive) external bus numbers to consecutive internal bus
%   numbers which start at 1. Changes are made to BUS, GEN, BRANCH and
%   optionally AREAS matrices, which are returned along with a vector of
%   indices I2E that can be passed to INT2EXT to perform the reverse
%   conversion, where EXTERNAL_BUS_NUMBER = I2E(INTERNAL_BUS_NUMBER)
%
%   Examples:
%       [i2e, bus, gen, branch, areas] = ext2int(bus, gen, branch, areas);
%       [i2e, bus, gen, branch] = ext2int(bus, gen, branch);
%
%   2.  MPC = EXT2INT(MPC)
%
%   If the input is a single MATPOWER case struct, then all isolated
%   buses, off-line generators and branches are removed along with any
%   generators, branches or areas connected to isolated buses. Then the
%   buses are renumbered consecutively, beginning at 1, and the
%   generators are sorted by increasing bus number. All of the related
%   indexing information and the original data matrices are stored in
%   an 'order' field in the struct to be used by INT2EXT to perform
%   the reverse conversions. If the case is already using internal
%   numbering it is returned unchanged.
%
%   Example:
%       mpc = ext2int(mpc);
%
%   3.  VAL = EXT2INT(MPC, VAL, ORDERING)
%       VAL = EXT2INT(MPC, VAL, ORDERING, DIM)
%       MPC = EXT2INT(MPC, FIELD, ORDERING)
%       MPC = EXT2INT(MPC, FIELD, ORDERING, DIM)
%
%   When given a case struct that has already been converted to
%   internal indexing, this function can be used to convert other data
%   structures as well by passing in 2 or 3 extra parameters in
%   addition to the case struct. If the value passed in the 2nd
%   argument is a column vector, it will be converted according to the
%   ORDERING specified by the 3rd argument (described below). If VAL
%   is an n-dimensional matrix, then the optional 4th argument (DIM,
%   default = 1) can be used to specify which dimension to reorder.
%   The return value in this case is the value passed in, converted
%   to internal indexing.
%
%   If the 2nd argument is a string or cell array of strings, it
%   specifies a field in the case struct whose value should be
%   converted as described above. In this case, the converted value
%   is stored back in the specified field, the original value is
%   saved for later use and the updated case struct is returned.
%   If FIELD is a cell array of strings, they specify nested fields.
%
%   The 3rd argument, ORDERING, is used to indicate whether the data
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
%       A_int = ext2int(mpc, A_ext, {'bus','bus','gen','gen'}, 2);
%
%       Converts an A matrix for user-supplied OPF constraints from
%       external to internal ordering, where the columns of the A
%       matrix correspond to bus voltage angles, then voltage
%       magnitudes, then generator real power injections and finally
%       generator reactive power injections.
%
%       gencost_int = ext2int(mpc, gencost_ext, {'gen','gen'}, 1);
%
%       Converts a GENCOST matrix that has both real and reactive power
%       costs (in rows 1--ng and ng+1--2*ng, respectively).
%
%       mpc = ext2int(mpc, {'reserves', 'cost'}, 'gen');
%
%       Reorders rows of mpc.reserves.cost to match internal generator
%       ordering.
%
%       mpc = ext2int(mpc, {'reserves', 'zones'}, 'gen', 2);
%
%       Reorders columns of mpc.reserves.zones to match internal
%       generator ordering.
%
%   The 'order' field of MPC used to store the indexing information
%   needed for subsequent internal to external conversion is structured
%   as:
%
%       order
%           state       'i' | 'e'
%           ext | int
%               areas
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
%           areas
%               status
%                   on
%                   off
%
%   See also INT2EXT.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2010 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if isstruct(bus)
    mpc = bus;
    if nargin == 1
        first = ~isfield(mpc, 'order');
        if first || mpc.order.state == 'e'
            %% define names for columns to data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE] = idx_bus;
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS] = idx_gen;
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS] = idx_brch;
            [AREA_I, PRICE_REF_BUS] = idx_area;

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
            if isfield(mpc, 'areas')
                if isempty(mpc.areas)           %% if areas field is empty
                    mpc = rmfield(mpc, 'areas');    %% delete it (so it gets ignored)
                else                            %% otherwise
                    o.ext.areas = mpc.areas;        %% save it
                end
            end

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
            if isfield(mpc, 'areas')
                as = bs(n2i(mpc.areas(:, PRICE_REF_BUS)));
                o.areas.status.on   = find(  as );
                o.areas.status.off  = find( ~as );
            end

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
            if isfield(mpc, 'areas') && ~isempty(o.areas.status.off)
                mpc.areas(o.areas.status.off, :) = [];
            end

            %% update size
            nb = size(mpc.bus, 1);

            %% apply consecutive bus numbering
            o.bus.i2e = mpc.bus(:, BUS_I);
            o.bus.e2i = sparse(max(o.bus.i2e), 1);
            o.bus.e2i(o.bus.i2e) = (1:nb)';
            mpc.bus(:, BUS_I)       = o.bus.e2i( mpc.bus(:, BUS_I)      );
            mpc.gen(:, GEN_BUS)     = o.bus.e2i( mpc.gen(:, GEN_BUS)    );
            mpc.branch(:, F_BUS)    = o.bus.e2i( mpc.branch(:, F_BUS)   );
            mpc.branch(:, T_BUS)    = o.bus.e2i( mpc.branch(:, T_BUS)   );
            if isfield(mpc, 'areas')
                mpc.areas(:, PRICE_REF_BUS) = o.bus.e2i( mpc.areas(:, PRICE_REF_BUS)  );
            end

            %% reorder gens in order of increasing bus number
            [tmp, o.gen.e2i] = sort(mpc.gen(:, GEN_BUS));
            [tmp, o.gen.i2e] = sort(o.gen.e2i);
            mpc.gen = mpc.gen(o.gen.e2i, :);

            if isfield(o, 'int')
                o = rmfield(o, 'int');
            end
            o.state = 'i';
            mpc.order = o;

            %% update gencost, A and N
            if isfield(mpc, 'gencost')
                ordering = {'gen'};         %% Pg cost only
                if size(mpc.gencost, 1) == 2*ng0
                    ordering{2} = 'gen';    %% include Qg cost
                end
                mpc = ext2int(mpc, 'gencost', ordering);
            end
            if isfield(mpc, 'A') || isfield(mpc, 'N')
                if dc
                    ordering = {'bus', 'gen'};
                else
                    ordering = {'bus', 'bus', 'gen', 'gen'};
                end
            end
            if isfield(mpc, 'A')
                mpc = ext2int(mpc, 'A', ordering, 2);
            end
            if isfield(mpc, 'N')
                mpc = ext2int(mpc, 'N', ordering, 2);
            end

            %% execute userfcn callbacks for 'ext2int' stage
            if isfield(mpc, 'userfcn')
                mpc = run_userfcn(mpc.userfcn, 'ext2int', mpc);
            end
        end

        i2e = mpc;
    else                    %% convert extra data
        ordering = branch;              %% rename argument
        if nargin < 4
            dim = 1;
        else
            dim = areas;                %% rename argument
        end
        if ischar(gen) || iscell(gen)   %% field
            field = gen;                %% rename argument
            if ischar(field)
                mpc.order.ext.(field) = mpc.(field);
                mpc.(field) = ext2int(mpc, mpc.(field), ordering, dim);
            else
                for k = 1:length(field)
                    s(k).type = '.';
                    s(k).subs = field{k};
                end
                mpc.order.ext = subsasgn(mpc.order.ext, s, subsref(mpc, s));
                mpc = subsasgn(mpc, s, ...
                    ext2int(mpc, subsref(mpc, s), ordering, dim));
            end
            i2e = mpc;
        else                            %% value
            val = gen;                  %% rename argument
            o = mpc.order;
            if ischar(ordering)         %% single set
                if strcmp(ordering, 'gen')
                    idx = o.(ordering).status.on(o.(ordering).e2i);
                else
                    idx = o.(ordering).status.on;
                end
                i2e = get_reorder(val, idx, dim);
            else                            %% multiple sets
                b = 0;  %% base
                for k = 1:length(ordering)
                    n = size(o.ext.(ordering{k}), 1);
                    v = get_reorder(val, b+(1:n), dim);
                    new_v{k} = ext2int(mpc, v, ordering{k}, dim);
                    b = b + n;
                end
                n = size(val, dim);
                if n > b                %% the rest
                    v = get_reorder(val, b+1:n, dim);
                    new_v{length(new_v)+1} = v;
                end
                i2e = cat(dim, new_v{:});
            end
        end
    end
else            %% old form
    %% define names for columns to data matrices
    [PQ, PV, REF, NONE, BUS_I] = idx_bus;
    [GEN_BUS] = idx_gen;
    [F_BUS, T_BUS] = idx_brch;
    [AREA_I, PRICE_REF_BUS] = idx_area;

    i2e = bus(:, BUS_I);
    e2i = sparse(max(i2e), 1);
    e2i(i2e) = (1:size(bus, 1))';

    bus(:, BUS_I)               = e2i( bus(:, BUS_I)            );
    gen(:, GEN_BUS)             = e2i( gen(:, GEN_BUS)          );
    branch(:, F_BUS)            = e2i( branch(:, F_BUS)         );
    branch(:, T_BUS)            = e2i( branch(:, T_BUS)         );
    if nargin > 3 && nargout > 4 && ~isempty(areas)
        areas(:, PRICE_REF_BUS) = e2i( areas(:, PRICE_REF_BUS)  );
    end
end
