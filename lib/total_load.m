function [Pd, Qd] = total_load(bus, gen, load_zone, opt)
%TOTAL_LOAD Returns vector of total load in each load zone.
%   PD = TOTAL_LOAD(BUS) returns active power demand for each zone.
%   PD = TOTAL_LOAD(BUS, GEN, LOAD_ZONE, OPT)
%   [PD, QD] = TOTAL_LOAD(...) returns both active and reative power
%   demand for each zone.
%
%   BUS - standard BUS matrix with nb rows, where the fixed active
%       and reactive loads are specified in columns PD and QD
%
%   GEN - (optional) standard GEN matrix with ng rows, where the
%       dispatchable loads are specified by columns PG, QG, PMIN,
%       QMIN and QMAX (in rows for which ISLOAD(GEN) returns true).
%       If GEN is empty, it assumes there are no dispatchable loads.
%
%   LOAD_ZONE - (optional) nb element vector where the value of
%       each element is either zero or the index of the load zone
%       to which the corresponding bus belongs. If LOAD_ZONE(b) = k
%       then the loads at bus b will added to the values of PD(k) and
%       QD(k). If LOAD_ZONE is empty, the default is defined as the areas
%       specified in the BUS matrix, i.e. LOAD_ZONE = BUS(:, BUS_AREA)
%       and load will have dimension = MAX(BUS(:, BUS_AREA)). If
%       LOAD_ZONE = 'all', the result is a scalar with the total system
%       load.
%
%   OPT - (optional) option struct, with the following fields:
%           'type'  -  string specifying types of loads to include, default
%                      is 'BOTH' if GEN is provided, otherwise 'FIXED'
%               'FIXED'        : sum only fixed loads
%               'DISPATCHABLE' : sum only dispatchable loads
%               'BOTH'         : sum both fixed and dispatchable loads
%           'nominal' -  1 : use nominal load for dispatchable loads
%                        0 : (default) use actual realized load for
%                             dispatchable loads
%
%       For backward compatibility with MATPOWER 4.x, OPT can also
%       take the form of a string, with the same options as OPT.type above.
%       In this case, again for backward compatibility, it is the "nominal"
%       load that is computed for dispatchable loads, not the actual
%       realized load. Using a string for OPT is deprecated and
%       will be removed in a future version.
%
%   Examples:
%       Return the total active load for each area as defined in BUS_AREA.
%
%       Pd = total_load(bus);
%
%       Return total active and reactive load, fixed and dispatchable, for
%       entire system.
%
%       [Pd, Qd] = total_load(bus, gen, 'all');
%
%       Return the total of the nominal dispatchable loads at buses 10-20.
%
%       load_zone = zeros(nb, 1);
%       load_zone(10:20) = 1;
%       opt = struct('type', 'DISPATCHABLE', 'nominal', 1);
%       Pd = total_load(bus, gen, load_zone, opt)
%
%   See also SCALE_LOAD.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2004-2014 by Power System Engineering Research Center (PSERC)
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

%% define constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
%% purposely being backward compatible with older MATPOWER
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, ...
    PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;

nb = size(bus, 1);          %% number of buses

%%-----  process inputs  -----
if nargin < 4
    opt = [];
    if nargin < 3
        load_zone = [];
        if nargin < 2
            gen = [];
        end
    end
end

%% default options
if ischar(opt)      %% convert old WHICH_TYPE string option to struct
    opt = struct('type', opt, 'nominal', 1);
else
    if ~isfield(opt, 'type')
        if isempty(gen)
            opt.type = 'FIXED';
        else
            opt.type = 'BOTH';
        end
    end
    if ~isfield(opt, 'nominal')
        opt.nominal = 0;
    end
end
switch upper(opt.type(1))
    case {'F', 'D', 'B'}
        %% OK
    otherwise
        error('total_load: OPT.type should be ''FIXED'', ''DISPATCHABLE'' or ''BOTH''');
end
want_Q      = (nargout > 1);
want_fixed  = (opt.type(1) == 'B' || opt.type(1) == 'F');
want_disp   = (opt.type(1) == 'B' || opt.type(1) == 'D');

%% initialize load_zone
if ischar(load_zone) && strcmp(load_zone, 'all')
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

    %% create map of external bus numbers to bus indices
    i2e = bus(:, BUS_I);
    e2i = sparse(max(i2e), 1);
    e2i(i2e) = (1:nb)';

    Cld = sparse(e2i(gen(:, GEN_BUS)), (1:ng)', is_ld, nb, ng);
    if opt.nominal      %% use nominal power
        Pdd = -Cld * gen(:, PMIN);      %% real power
        if want_Q
            Q = zeros(ng, 1);
            Q(ld) = (gen(ld, QMIN) == 0) .* gen(ld, QMAX) + ...
                    (gen(ld, QMAX) == 0) .* gen(ld, QMIN);
            Qdd = -Cld * Q;             %% reactive power
        end
    else                %% use realized actual power dispatch
        Pdd = -Cld * gen(:, PG);        %% real power
        if want_Q
            Qdd = -Cld * gen(:, QG);    %% reactive power
        end
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
    idx = find( load_zone == k & bus(:, BUS_TYPE) ~= NONE);
    Pd(k) = sum(Pdf(idx)) + sum(Pdd(idx));
    if want_Q
        Qd(k) = sum(Qdf(idx)) + sum(Qdd(idx));
    end
end
