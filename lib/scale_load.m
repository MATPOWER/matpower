function [bus, gen, gencost] = scale_load(dmd, bus, gen, load_zone, opt, gencost)
%SCALE_LOAD Scales fixed and/or dispatchable loads.
%   MPC = SCALE_LOAD(LOAD, MPC);
%   MPC = SCALE_LOAD(LOAD, MPC, LOAD_ZONE)
%   MPC = SCALE_LOAD(LOAD, MPC, LOAD_ZONE, OPT)
%   BUS = SCALE_LOAD(LOAD, BUS);
%   [BUS, GEN] = SCALE_LOAD(LOAD, BUS, GEN, LOAD_ZONE, OPT)
%   [BUS, GEN, GENCOST] = ...
%       SCALE_LOAD(LOAD, BUS, GEN, LOAD_ZONE, OPT, GENCOST)
%
%   Scales active (and optionally reactive) loads in each zone by a
%   zone-specific ratio, i.e. R(k) for zone k. Inputs are ...
%
%   LOAD - Each element specifies the amount of scaling for the
%       corresponding load zone, either as a direct scale factor
%       or as a target quantity. If there are nz load zones this
%       vector has nz elements.
%
%   MPC - standard MATPOWER case struct or case file name
%       
%   BUS - standard BUS matrix with nb rows, where the fixed active
%       and reactive loads available for scaling are specified in
%       columns PD and QD
%
%   GEN - (optional) standard GEN matrix with ng rows, where the
%       dispatchable loads available for scaling are specified by
%       columns PG, QG, PMIN, QMIN and QMAX (in rows for which
%       ISLOAD(GEN) returns true). If GEN is empty, it assumes
%       there are no dispatchable loads.
%
%   LOAD_ZONE - (optional) nb element vector where the value of
%       each element is either zero or the index of the load zone
%       to which the corresponding bus belongs. If LOAD_ZONE(b) = k
%       then the loads at bus b will be scaled according to the
%       value of LOAD(k). If LOAD_ZONE(b) = 0, the loads at bus b
%       will not be modified. If LOAD_ZONE is empty, the default is
%       determined by the dimensions of the LOAD vector. If LOAD is
%       a scalar, a single system-wide zone including all buses is
%       used, i.e. LOAD_ZONE = ONES(nb, 1). If LOAD is a vector, the
%       default LOAD_ZONE is defined as the areas specified in the
%       BUS matrix, i.e. LOAD_ZONE = BUS(:, BUS_AREA), and LOAD
%       should have dimension = MAX(BUS(:, BUS_AREA)).
%
%   OPT - (optional) struct with three possible fields, 'scale',
%       'pq' and 'which' that determine the behavior as follows:
%
%     OPT.scale (default is 'FACTOR')
%       'FACTOR'   : LOAD consists of direct scale factors, where
%                    LOAD(k) = scale factor R(k) for zone k
%       'QUANTITY' : LOAD consists of target quantities, where
%                    LOAD(k) = desired total active load in MW for
%                    zone k after scaling by an appropriate R(k)
%
%     OPT.pq    (default is 'PQ')
%       'PQ' : scale both active and reactive loads
%       'P'  : scale only active loads
%
%     OPT.which (default is 'BOTH' if GEN is provided, else 'FIXED')
%       'FIXED'        : scale only fixed loads
%       'DISPATCHABLE' : scale only dispatchable loads
%       'BOTH'         : scale both fixed and dispatchable loads
%
%     OPT.cost : (default = -1) flag to include cost in scaling or not
%       -1 : include cost if gencost is available
%        0 : do not include cost
%        1 : include cost (error if gencost not available)
%
%   GENCOST - (optional) standard GENCOST matrix with ng (or 2*ng)
%       rows, where the dispatchable load rows are determined by
%       the GEN matrix. If included, the quantity axis of the marginal
%       "cost" or benefit function of any dispatchable loads will be
%       scaled with the size of the load itself (using MODCOST twice,
%       once with MODTYPE equal to SCALE_F and once with SCALE_X).
%
%   Examples:
%       Scale all real and reactive fixed loads up by 10%.
%
%       bus = scale_load(1.1, bus);
%
%       Scale all active loads (fixed and dispatchable) at the first 10
%       buses so their total equals 100 MW, and at next 10 buses so their
%       total equals 50 MW.
%
%       load_zone = zeros(nb, 1);
%       load_zone(1:10) = 1;
%       load_zone(11:20) = 2;
%       opt = struct('pq', 'P', 'scale', 'QUANTITY');
%       dmd = [100; 50];
%       [bus, gen] = scale_load(dmd, bus, gen, load_zone, opt);
%
%   See also TOTAL_LOAD.

%   MATPOWER
%   Copyright (c) 2004-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
%% purposely being backward compatible with older MATPOWER
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, ...
    PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;

%%-----  process inputs  -----
if ischar(bus)      %% passing in case as file name string
    bus = loadcase(bus);
end
if isstruct(bus)
    use_mpc = 1;
    if nargin < 4
        load_zone = struct;
        if nargin < 3
            gen = [];
        end
    end
    %% shift and reassign inputs
    opt = load_zone;
    load_zone = gen;
    mpc = bus;
    gen = mpc.gen;
    bus = mpc.bus;
    if isfield(mpc, 'gencost')
        gencost = mpc.gencost;
    else
        gencost = [];
    end
    if nargout > 1
        error('scale_load: too many output arguments')
    end
else
    use_mpc = 0;
    if nargin < 6
        gencost = [];
        if nargin < 5
            opt = struct;
            if nargin < 4
                load_zone = [];
                if nargin < 3
                    gen = [];
                end
            end
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
if ~isfield(opt, 'cost')
    opt.cost = -1;          %% -1, 0, or 1
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
if opt.cost == -1
    if isempty(gencost)
        opt.cost = 0;
    else
        opt.cost = 1;
    end
end
if ~use_mpc && nargout < 3 && opt.cost
    error('scale_load: missing gencost as output argument');
end
if nargout > 2 && isempty(gencost)
    error('scale_load: missing gencost as input argument');
end

%% create dispatchable load connection matrix
nb = size(bus, 1);          %% number of buses
if ~isempty(gen)
    ng = size(gen, 1);
    is_ld = isload(gen) & gen(:, GEN_STATUS) > 0;
    ld = find(is_ld);

    %% create map of external bus numbers to bus indices
    i2e = bus(:, BUS_I);
    e2i = sparse(max(i2e), 1);
    e2i(i2e) = (1:nb)';

    Cld = sparse(e2i(gen(:, GEN_BUS)), (1:ng)', is_ld, nb, ng);
else
    ng = [];
    ld = [];
end

if isempty(load_zone)
    if length(dmd) == 1         %% make a single zone of all load buses
        load_zone = zeros(nb, 1);                           %% initialize
        load_zone(bus(:, PD) ~= 0 | bus(:, QD) ~= 0) = 1;   %% FIXED loads
        if ~isempty(gen)
            load_zone(e2i(gen(ld, GEN_BUS))) = 1;   %% DISPATCHABLE loads
        end
    else                        %% use areas defined in bus data as zones
        load_zone = bus(:, BUS_AREA);
    end
end

%% check load_zone to make sure it's consistent with size of load vector
if max(load_zone) > length(dmd)
    error('scale_load: load vector must have a value for each load zone specified');
end

%%-----  compute scale factors for each zone  -----
scale = dmd;
Pdd = zeros(nb, 1);     %% dispatchable P at each bus
if opt.scale(1) == 'Q'  %% 'QUANTITY'
    %% find load capacity from dispatchable loads
    if ~isempty(gen)
        Pdd = -Cld * gen(:, PMIN);
    end

    %% compute scale factors
    for k = 1:length(dmd)
        idx = find( load_zone == k );
        fixed = sum(bus(idx, PD));
        dispatchable = sum(Pdd(idx));
        total = fixed + dispatchable;
        if opt.which(1) == 'B'      %% 'BOTH'
            if total ~= 0
                scale(k) = dmd(k) / total;
            elseif dmd(k) == total
                scale(k) = 1;
            else
                error('scale_load: impossible to make zone %d load equal %g by scaling non-existent loads', k, dmd(k));
            end
        elseif opt.which(1) == 'F'  %% 'FIXED'
            if fixed ~= 0
                scale(k) = (dmd(k) - dispatchable) / fixed;
            elseif dmd(k) == dispatchable
                scale(k) = 1;
            else
                error('scale_load: impossible to make zone %d load equal %g by scaling non-existent fixed load', k, dmd(k));
            end
        elseif opt.which(1) == 'D'  %% 'DISPATCHABLE'
            if dispatchable ~= 0
                scale(k) = (dmd(k) - fixed) / dispatchable;
            elseif dmd(k) == fixed
                scale(k) = 1;
            else
                error('scale_load: impossible to make zone %d load equal %g by scaling non-existent dispatchable load', k, dmd(k));
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
        if strcmp(opt.pq, 'PQ')
            bus(idx, QD) = bus(idx, QD) * scale(k);
        end
    end
end

%% dispatchable loads
if opt.which(1) ~= 'F'      %% includes 'DISPATCHABLE', not 'FIXED' only
    for k = 1:length(scale)
        idx = find( load_zone == k );
        [junk, i, junk2] = intersect(e2i(gen(ld, GEN_BUS)), idx);
        ig = ld(i);

        gen(ig, [PG PMIN]) = gen(ig, [PG PMIN]) * scale(k);
        if opt.cost
            gencost(ig, :) = modcost(gencost(ig, :), scale(k), 'SCALE_F');
            gencost(ig, :) = modcost(gencost(ig, :), scale(k), 'SCALE_X');
        end
        if strcmp(opt.pq, 'PQ')
            gen(ig, [QG QMIN QMAX]) = gen(ig, [QG QMIN QMAX]) * scale(k);
            if opt.cost
                [pcost, qcost] = pqcost(gencost, ng);
                if ~isempty(qcost)
                    qcost(ig, :) = modcost(qcost(ig, :), scale(k), 'SCALE_F');
                    qcost(ig, :) = modcost(qcost(ig, :), scale(k), 'SCALE_X');
                    gencost = [pcost; qcost];
                end
            end
        end
    end
end

%% re-package outputs if necessary
if use_mpc
    mpc.bus = bus;
    mpc.gen = gen;
    if opt.cost
        mpc.gencost = gencost;
    end
    bus = mpc;
end
