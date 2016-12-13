function [ref, pv, pq] = bustypes(bus, gen)
%BUSTYPES   Builds index lists for each type of bus (REF, PV, PQ).
%   [REF, PV, PQ] = BUSTYPES(BUS, GEN)
%   Generators with "out-of-service" status are treated as PQ buses with
%   zero generation (regardless of Pg/Qg values in gen). Expects BUS and
%   GEN have been converted to use internal consecutive bus numbering.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% get generator status
% bus_gen_status = zeros(size(bus, 1), 1);
% bus_gen_status(gen(:, GEN_BUS)) = gen(:, GEN_STATUS) > 0;
nb = size(bus, 1);
ng = size(gen, 1);
Cg = sparse(gen(:, GEN_BUS), (1:ng)', gen(:, GEN_STATUS) > 0, nb, ng);  %% gen connection matrix
                                        %% element i, j is 1 if, generator j at bus i is ON
bus_gen_status = Cg * ones(ng, 1);      %% number of generators at each bus that are ON


%% form index lists for slack, PV, and PQ buses
ref = find(bus(:, BUS_TYPE) == REF & bus_gen_status);   %% reference bus index
pv  = find(bus(:, BUS_TYPE) == PV  & bus_gen_status);   %% PV bus indices
pq  = find(bus(:, BUS_TYPE) == PQ | ~bus_gen_status);   %% PQ bus indices

%% pick a new reference bus if for some reason there is none (may have been shut down)
if isempty(ref)
%     if isempty(pv)
%         %% no PV bus left to convert to reference bus
%     else
        ref = pv(1);    %% use the first PV bus
        pv(1) = [];     %% delete it from PV list
%     end
end
