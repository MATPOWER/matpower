function [Pd, Qd] = total_load(bus, gen, load_zone, which_type)
%TOTAL_LOAD Returns vector of total load in each load zone.
%
%   Pd = total_load(bus)
%   Pd = total_load(bus, gen)
%   Pd = total_load(bus, gen, load_zone)
%   Pd = total_load(bus, gen, load_zone, which_type)
%
%   [Pd, Qd] = total_load(bus)
%   [Pd, Qd] = total_load(bus, gen)
%   [Pd, Qd] = total_load(bus, gen, load_zone)
%   [Pd, Qd] = total_load(bus, gen, load_zone, which_type)
%
%
%   load_zone - (optional) nb element vector where the value of
%       each element is either zero or the index of the load zone
%       to which the corresponding bus belongs. If load_zone(b) = k
%       then the loads at bus b will added to the value of load(k).
%       If load_zone is empty, the default is defined as the areas
%       specified in the bus matrix, i.e. load_zone = bus(:, BUS_AREA)
%       and load will have dimension = max(bus(:, BUS_AREA)). If
%       load_zone = 'all', the result is a scalar with the total system
%       load.
%
%   which_type - (optional) (default is 'BOTH' if gen is provided, else 'FIXED')
%       'FIXED'        : sum only fixed loads
%       'DISPATCHABLE' : sum only dispatchable loads
%       'BOTH'         : sum both fixed and dispatchable loads
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
if nargin < 4
    which_type = [];
    if nargin < 3
        load_zone = [];
        if nargin < 2
            gen = [];
        end
    end
end

%% fill out and check which_type
if isempty(gen)
    which_type = 'FIXED';
end
if isempty(which_type) && ~isempty(gen)
    which_type = 'BOTH';     %% 'FIXED', 'DISPATCHABLE' or 'BOTH'
end
if which_type(1) ~= 'F' && which_type(1) ~= 'D' && which_type(1) ~= 'B'
    error('total_load: which_type should be ''FIXED'', ''DISPATCHABLE'' or ''BOTH''');
end
want_Q      = (nargout > 1);
want_fixed  = (which_type(1) == 'B' || which_type(1) == 'F');
want_disp   = (which_type(1) == 'B' || which_type(1) == 'D');

%% initialize load_zone
if isstr(load_zone) && strcmp(load_zone, 'all')
    load_zone = ones(nb, 1);        %% make a single zone of all buses
elseif isempty(load_zone)
    load_zone = bus(:, BUS_AREA);   %% use areas defined in bus data as zones
end
nz = max(load_zone);    %% number of load zones

%% fixed load at each bus, & initialize dispatchable
if want_fixed
    Pdf = bus(:, PD);       %% real power
    if want_Q
        Qdf = bus(:, QD);   %% reactive power
    end
else
    Pdf = zeros(nb, 1);     %% real power
    if want_Q
        Qdf = zeros(nb, 1); %% reactive power
    end
end

%% dispatchable load at each bus 
if want_disp            %% need dispatchable
    ng = size(gen, 1);
    is_ld = isload(gen) & gen(:, GEN_STATUS) > 0;
    ld = find(is_ld);
    Cld = sparse(gen(:, GEN_BUS), [1:ng]', is_ld, nb, ng);
    Pdd = -Cld * gen(:, PMIN);      %% real power
    if want_Q
        Q = zeros(ng, 1);
        Q(ld) = (gen(ld, QMIN) == 0) .* gen(ld, QMAX) + ...
                (gen(ld, QMAX) == 0) .* gen(ld, QMIN);
        Qdd = -Cld * Q;             %% reactive power
    end
else
    Pdd = zeros(nb, 1);
    if want_Q
        Qdd = zeros(nb, 1);
    end
end

%% compute load sums
Pd = zeros(nz, 1);
if want_Q
    Qd = zeros(nz, 1);
end
for k = 1:nz
    idx = find( load_zone == k );
    Pd(k) = sum(Pdf(idx)) + sum(Pdd(idx));
    if want_Q
        Qd(k) = sum(Qdf(idx)) + sum(Qdd(idx));
    end
end

return;
