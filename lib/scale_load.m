function [bus, gen] = scale_load(load, bus, gen, load_zone, opt)
%SCALE_LOAD Scales fixed and/or dispatchable loads.
%
%   bus = scale_load(load, bus)
%   [bus, gen] = scale_load(load, bus, gen)
%   [bus, gen] = scale_load(load, bus, gen, load_zone)
%   [bus, gen] = scale_load(load, bus, gen, load_zone, opt)
%
%   Scales active (and optionally reactive) loads in each zone by a
%   zone-specific ratio, i.e. R(k) for zone k. Inputs are ...
%
%   load - Each element specifies the amount of scaling for the
%       corresponding load zone, either as a direct scale factor
%       or as a target quantity. If there are nz load zones this
%       vector has nz elements.
%
%   bus - standard bus matrix with nb rows, where the fixed active
%       and reactive loads available for scaling are specified in
%       columns PD and QD
%
%   gen - (optional) standard gen matrix with ng rows, where the
%       dispatchable loads available for scaling are specified by
%       columns PG, QG, PMIN, QMIN and QMAX (in rows for which
%       isload(gen) returns true). If gen is empty, it assumes
%       there are no dispatchable loads.
%
%   load_zone - (optional) nb element vector where the value of
%       each element is either zero or the index of the load zone
%       to which the corresponding bus belongs. If load_zone(b) = k
%       then the loads at bus b will be scaled according to the
%       value of load(k). If load_zone(b) = 0, the loads at bus b
%       will not be modified. If load_zone is empty, the default is
%       determined by the dimensions of the load vector. If load is
%       a scalar, a single system-wide zone including all buses is
%       used, i.e. loadzone = ones(nb, 1). If load is a vector, the
%       default load_zone is defined as the areas specified in the
%       bus matrix, i.e. load_zone = bus(:, BUS_AREA), and load
%       should have dimension = max(bus(:, BUS_AREA)).
%
%   opt - (optional) struct with three possible fields, 'scale',
%       'pq' and 'which' that determine the behavior as follows:
%
%     opt.scale (default is 'FACTOR')
%       'FACTOR'    : load consists of direct scale factors, where
%                    load(k) = scale factor R(k) for zone k
%       'QUANTITY' : load consists of target quantities, where
%                    load(k) = desired total active load for zone k
%                    after scaling by an appropriate R(k)
%
%     opt.pq    (default is 'PQ')
%       'PQ' : scale both active and reactive loads
%       'P'  : scale only active loads
%
%     opt.which (default is 'BOTH' if gen is provided, else 'FIXED')
%       'FIXED'        : scale only fixed loads
%       'DISPATCHABLE' : scale only dispatchable loads
%       'BOTH'         : scale both fixed and dispatchable loads
%
%   Assumes consecutive bus numbering when dealing with dispatchable loads.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2004-2009 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
%% purposely being backward compatible with older MATPOWER
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, ...
    PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch;

nb = size(bus, 1);          %% number of buses

%%-----  process inputs  -----
if nargin < 5
    opt = struct;
    if nargin < 4
        load_zone = [];
        if nargin < 3
            gen = [];
        end
    end
end

%% fill out and check opt
if isempty(gen)
    opt.which = 'FIXED';
end
if ~isfield(opt, 'pq')
    opt.pq = 'PQ';          %% 'PQ' or 'P'
end
if ~isfield(opt, 'which')
    opt.which = 'BOTH';     %% 'FIXED', 'DISPATCHABLE' or 'BOTH'
end
if ~isfield(opt, 'scale')
    opt.scale = 'FACTOR';   %% 'FACTOR' or 'QUANTITY'
end
if ~strcmp(opt.pq, 'P') && ~strcmp(opt.pq, 'PQ')
    error('scale_load: opt.pq must equal ''PQ'' or ''P''');
end
if opt.which(1) ~= 'F' && opt.which(1) ~= 'D' && opt.which(1) ~= 'B'
    error('scale_load: opt.which should be ''FIXED'', ''DISPATCHABLE'' or ''BOTH''');
end
if opt.scale(1) ~= 'F' && opt.scale(1) ~= 'Q'
    error('scale_load: opt.scale should be ''FACTOR'' or ''QUANTITY''');
end
if isempty(gen) && opt.which(1) ~= 'F'
    error('scale_load: need gen matrix to scale dispatchable loads');
end

%% create dispatchable load connection matrix
if ~isempty(gen)
    ng = size(gen, 1);
    is_ld = isload(gen) & gen(:, GEN_STATUS) > 0;
    ld = find(is_ld);
    Cld = sparse(gen(:, GEN_BUS), (1:ng)', is_ld, nb, ng);
else
    ng = [];
    ld = [];
end

if isempty(load_zone)
    if length(load) == 1        %% make a single zone of all load buses
        load_zone = zeros(nb, 1);   %% initialize
        k = find(bus(:, PD));
        load_zone(k) = 1;                       %% FIXED loads
        if ~isempty(gen)
            load_zone(gen(ld, GEN_BUS)) = 1;    %% DISPATCHABLE loads
        end
    else                        %% use areas defined in bus data as zones
        load_zone = bus(:, BUS_AREA);
    end
end

%% check load_zone to make sure it's consistent with size of load vector
if max(load_zone) > length(load)
    error('scale_load: load vector must have a value for each load zone specified');
end

%%-----  compute scale factors for each zone  -----
scale = load;
Pdd = zeros(nb, 1);     %% dispatchable P at each bus
Qdd = zeros(nb, 1);     %% dispatchable Q at each bus
if opt.scale(1) == 'Q'  %% 'QUANTITY'
    %% find load capacity from dispatchable loads
    if ~isempty(gen)
        Q = zeros(ng, 1);
        Q(ld) = (gen(ld, QMIN) == 0) .* gen(ld, QMAX) + ...
                (gen(ld, QMAX) == 0) .* gen(ld, QMIN);
        Pdd = -Cld * gen(:, PMIN);
        Qdd = -Cld * Q;
    end

    %% compute scale factors
    for k = 1:length(load)
        idx = find( load_zone == k );
        fixed = sum(bus(idx, PD));
        dispatchable = sum(Pdd(idx));
        total = fixed + dispatchable;
        if opt.which(1) == 'B'      %% 'BOTH'
            if total ~= 0
                scale(k) = load(k) / total;
            elseif load(k) == total
                scale(k) = 1;
            else
                error('scale_load: impossible to make zone %d load equal %g by scaling non-existent loads', k, load(k));
            end
        elseif opt.which(1) == 'F'  %% 'FIXED'
            if fixed ~= 0
                scale(k) = (load(k) - dispatchable) / fixed;
            elseif load(k) == dispatchable
                scale(k) = 1;
            else
                error('scale_load: impossible to make zone %d load equal %g by scaling non-existent fixed load', k, load(k));
            end
        elseif opt.which(1) == 'D'  %% 'DISPATCHABLE'
            if dispatchable ~= 0
                scale(k) = (load(k) - fixed) / dispatchable;
            elseif load(k) == fixed
                scale(k) = 1;
            else
                error('scale_load: impossible to make zone %d load equal %g by scaling non-existent dispatchable load', k, load(k));
            end
        end
    end
end

%%-----  do the scaling  -----
%% fixed loads
if opt.which(1) ~= 'D'      %% includes 'FIXED', not 'DISPATCHABLE' only
    for k = 1:length(scale)
        idx = find( load_zone == k );
        bus(idx, PD) = bus(idx, PD) * scale(k);
        if opt.pq == 'PQ'
            bus(idx, QD) = bus(idx, QD) * scale(k);
        end
    end
end

%% dispatchable loads
if opt.which(1) ~= 'F'      %% includes 'DISPATCHABLE', not 'FIXED' only
    for k = 1:length(scale)
        idx = find( load_zone == k );
        [junk, i, junk2] = intersect(gen(ld, GEN_BUS), idx);
        ig = ld(i);

        gen(ig, [PG PMIN]) = gen(ig, [PG PMIN]) * scale(k);
        if opt.pq == 'PQ'
            gen(ig, [QG QMIN QMAX]) = gen(ig, [QG QMIN QMAX]) * scale(k);
        end
    end
end
