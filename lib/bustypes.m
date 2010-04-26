function [ref, pv, pq] = bustypes(bus, gen)
%BUSTYPES   Builds index lists for each type of bus (REF, PV, PQ).
%   [REF, PV, PQ] = BUSTYPES(BUS, GEN)
%   Generators with "out-of-service" status are treated as PQ buses with
%   zero generation (regardless of Pg/Qg values in gen). Expects BUS and
%   GEN have been converted to use internal consecutive bus numbering.

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
    ref = pv(1);                %% use the first PV bus
    pv = pv(2:length(pv));      %% take it off PV list
end
