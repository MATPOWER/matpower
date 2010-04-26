function t_scale_load(quiet)
%T_SCALE_LOAD  Tests for code in SCALE_LOAD.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008-2010 by Power System Engineering Research Center (PSERC)
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

if nargin < 1
    quiet = 0;
end

n_tests = 275;

t_begin(n_tests, quiet);

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

mpc = loadcase('t_auction_case');
mpc.gen(8, GEN_BUS) = 2;    %% multiple d. loads per area, same bus as gen
mpc.gen(8, [QG QMIN QMAX]) = [ 3 0 3 ];
%% put it load before gen in matrix
mpc.gen = [mpc.gen(8, :); mpc.gen(1:7, :); mpc.gen(9, :)];
ld = find(isload(mpc.gen));
for k = 1:3
    a{k} = find(mpc.bus(:, BUS_AREA) == k); %% buses in area k
    [junk, tmp, junk2] = intersect(mpc.gen(ld, GEN_BUS), a{k});
    lda{k} = ld(tmp);                       %% disp loads in area k
end
for k = 1:3
    area(k).fixed.p = sum(mpc.bus(a{k}, PD));
    area(k).fixed.q = sum(mpc.bus(a{k}, QD));
    area(k).disp.p = -sum(mpc.gen(lda{k}, PMIN));
    area(k).disp.qmin = -sum(mpc.gen(lda{k}, QMIN));
    area(k).disp.qmax = -sum(mpc.gen(lda{k}, QMAX));
    area(k).disp.q = area(k).disp.qmin + area(k).disp.qmax;
    area(k).both.p = area(k).fixed.p + area(k).disp.p;
    area(k).both.q = area(k).fixed.q + area(k).disp.q;
end
total.fixed.p = sum(mpc.bus(:, PD));
total.fixed.q = sum(mpc.bus(:, QD));
total.disp.p = -sum(mpc.gen(ld, PMIN));
total.disp.qmin = -sum(mpc.gen(ld, QMIN));
total.disp.qmax = -sum(mpc.gen(ld, QMAX));
total.disp.q = total.disp.qmin + total.disp.qmax;
total.both.p = total.fixed.p + total.disp.p;
total.both.q = total.fixed.q + total.disp.q;

%%-----  single load zone, one scale factor  -----
load = 2;
t = 'all fixed loads (PQ) * 2 : ';
bus = scale_load(load, mpc.bus);
t_is(sum(bus(:, PD)), load*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), load*total.fixed.q, 8, [t 'total fixed Q']);
opt = struct('which', 'FIXED');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), load*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), load*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all fixed loads (P) * 2 : ';
opt = struct('pq', 'P');
bus = scale_load(load, mpc.bus, [], [], opt);
t_is(sum(bus(:, PD)), load*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
opt = struct('pq', 'P', 'which', 'FIXED');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), load*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all loads (PQ) * 2 : ';
[bus, gen] = scale_load(load, mpc.bus, mpc.gen);
t_is(sum(bus(:, PD)), load*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), load*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), load*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), load*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), load*total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all loads (P) * 2 : ';
opt = struct('pq', 'P');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), load*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), load*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all disp loads (PQ) * 2 : ';
opt = struct('which', 'DISPATCHABLE');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), load*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), load*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), load*total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all disp loads (P) * 2 : ';
opt = struct('pq', 'P', 'which', 'DISPATCHABLE');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), load*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);

%%-----  single load zone, one scale quantity  -----
load = 200;
t = 'all fixed loads (PQ) => total = 200 : ';
opt = struct('scale', 'QUANTITY');
bus = scale_load(load, mpc.bus, [], [], opt);
t_is(sum(bus(:, PD)), load, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), load/total.fixed.p*total.fixed.q, 8, [t 'total fixed Q']);
opt = struct('scale', 'QUANTITY', 'which', 'FIXED');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), load-total.disp.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), (load-total.disp.p)/total.fixed.p*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all fixed loads (P) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'pq', 'P');
bus = scale_load(load, mpc.bus, [], [], opt);
t_is(sum(bus(:, PD)), load, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
opt = struct('scale', 'QUANTITY', 'pq', 'P', 'which', 'FIXED');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), load-total.disp.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all loads (PQ) => total = 200 : ';
opt = struct('scale', 'QUANTITY');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), load/total.both.p*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), load/total.both.p*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), load/total.both.p*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), load/total.both.p*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), load/total.both.p*total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all loads (P) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'pq', 'P');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), load/total.both.p*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), load/total.both.p*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all disp loads (PQ) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'which', 'DISPATCHABLE');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), load-total.fixed.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), (load-total.fixed.p)/total.disp.p*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), (load-total.fixed.p)/total.disp.p*total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all disp loads (P) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'pq', 'P', 'which', 'DISPATCHABLE');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), load-total.fixed.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);

%%-----  3 zones, area scale factors  -----
t = 'area fixed loads (PQ) * [3 2 1] : ';
load = [3 2 1];
bus = scale_load(load, mpc.bus);
for k = 1:length(load)
    t_is(sum(bus(a{k}, PD)), load(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), load(k)*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
end
opt = struct('which', 'FIXED');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
for k = 1:length(load)
    t_is(sum(bus(a{k}, PD)), load(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), load(k)*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

t = 'area fixed loads (P) * [3 2 1] : ';
load = [3 2 1];
opt = struct('pq', 'P');
bus = scale_load(load, mpc.bus, [], [], opt);
for k = 1:length(load)
    t_is(sum(bus(a{k}, PD)), load(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
end
opt = struct('pq', 'P', 'which', 'FIXED');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
for k = 1:length(load)
    t_is(sum(bus(a{k}, PD)), load(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

t = 'all area loads (PQ) * [3 2 1] : ';
[bus, gen] = scale_load(load, mpc.bus, mpc.gen);
for k = 1:length(load)
    t_is(sum(bus(a{k}, PD)), load(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), load(k)*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), load(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), load(k)*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), load(k)*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

t = 'all area loads (P) * [3 2 1] : ';
opt = struct('pq', 'P');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
for k = 1:length(load)
    t_is(sum(bus(a{k}, PD)), load(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), load(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

t = 'area disp loads (PQ) * [3 2 1] : ';
opt = struct('which', 'DISPATCHABLE');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
for k = 1:length(load)
    t_is(sum(bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), load(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), load(k)*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), load(k)*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

t = 'area disp loads (P) * [3 2 1] : ';
opt = struct('pq', 'P', 'which', 'DISPATCHABLE');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
for k = 1:length(load)
    t_is(sum(bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), load(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

%%-----  3 zones, area scale quantities  -----
t = 'area fixed loads (PQ) => total = [100 80 60] : ';
load = [100 80 60];
opt = struct('scale', 'QUANTITY');
bus = scale_load(load, mpc.bus, [], [], opt);
for k = 1:length(load)
    t_is(sum(bus(a{k}, PD)), load(k), 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), load(k)/area(k).fixed.p*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
end
opt = struct('scale', 'QUANTITY', 'which', 'FIXED');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
for k = 1:length(load)
    t_is(sum(bus(a{k}, PD)), load(k)-area(k).disp.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), (load(k)-area(k).disp.p)/area(k).fixed.p*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

t = 'area fixed loads (P) => total = [100 80 60] : ';
load = [100 80 60];
opt = struct('scale', 'QUANTITY', 'pq', 'P');
bus = scale_load(load, mpc.bus, [], [], opt);
for k = 1:length(load)
    t_is(sum(bus(a{k}, PD)), load(k), 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
end
opt = struct('scale', 'QUANTITY', 'pq', 'P', 'which', 'FIXED');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
for k = 1:length(load)
    t_is(sum(bus(a{k}, PD)), load(k)-area(k).disp.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

t = 'all area loads (PQ) => total = [100 80 60] : ';
opt = struct('scale', 'QUANTITY');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
for k = 1:length(load)
    t_is(sum(bus(a{k}, PD)), load(k)/area(k).both.p*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), load(k)/area(k).both.p*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), load(k)/area(k).both.p*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), load(k)/area(k).both.p*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), load(k)/area(k).both.p*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

t = 'all area loads (P) => total = [100 80 60] : ';
opt = struct('scale', 'QUANTITY', 'pq', 'P');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
for k = 1:length(load)
    t_is(sum(bus(a{k}, PD)), load(k)/area(k).both.p*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), load(k)/area(k).both.p*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

t = 'area disp loads (PQ) => total = [100 80 60] : throws expected exception';
load = [100 80 60];
opt = struct('scale', 'QUANTITY', 'which', 'DISPATCHABLE');
err = 0;
try
    [bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
catch
    [msg, id] = lasterr;
    expected = 'scale_load: impossible to make zone 2 load equal 80 by scaling non-existent dispatchable load';
    if ~isempty(findstr(expected, msg))
        err = 1;
    end
end
t_ok(err, t);

t = 'area disp loads (PQ) => total = [100 74.3941 60] : ';
load = [100 area(2).fixed.p 60];
opt = struct('scale', 'QUANTITY', 'which', 'DISPATCHABLE');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
for k = 1:length(load)
    t_is(sum(bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), load(k)-area(k).fixed.p, 8, sprintf('%s area %d disp P', t, k));
    if k == 2
        t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    else
        t_is(-sum(gen(lda{k}, QMIN)), (load(k)-area(k).fixed.p)/area(k).disp.p*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(gen(lda{k}, QMAX)), (load(k)-area(k).fixed.p)/area(k).disp.p*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

t = 'area disp loads (P) => total = [100 74.3941 60] : ';
opt = struct('scale', 'QUANTITY', 'pq', 'P', 'which', 'DISPATCHABLE');
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, [], opt);
for k = 1:length(load)
    t_is(sum(bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), load(k)-area(k).fixed.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

%%-----  explict single load zone  -----
t = 'explicit single load zone';
load_zone = zeros(1, size(mpc.bus, 1));
load_zone([3 4]) = 1;
load = 2;
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, load_zone);
Pd = mpc.bus(:, PD);
Pd([3 4]) = load * Pd([3 4]);
t_is( bus(:, PD), Pd, 8, t);

%%-----  explict multiple load zone  -----
t = 'explicit multiple load zone';
load_zone = zeros(1, size(mpc.bus, 1));
load_zone([3 4]) = 1;
load_zone([7 8]) = 2;
load = [2 0.5];
[bus, gen] = scale_load(load, mpc.bus, mpc.gen, load_zone);
Pd = mpc.bus(:, PD);
Pd([3 4]) = load(1) * Pd([3 4]);
Pd([7 8]) = load(2) * Pd([7 8]);
t_is( bus(:, PD), Pd, 8, t);

t_end;
