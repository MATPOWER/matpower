function t_scale_load(quiet)
%T_SCALE_LOAD  Tests for code in SCALE_LOAD.

%   MATPOWER
%   Copyright (c) 2008-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

n_tests = 846;

t_begin(n_tests, quiet);

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

mpc = loadcase('t_auction_case');
mpc.gen(8, GEN_BUS) = 2;    %% multiple d. loads per area, same bus as gen
mpc.gen(8, [QG QMIN QMAX]) = [ 3 0 3 ];
mpc.gencost(7, COST:end) = [-30 -600 -20 -300 -10 -100 0 0];    % 10, 20, 30
mpc.gencost(8, COST:end) = [-30  -60 -20  -30 -10  -10 0 0];    % 1, 2, 3
mpc.gencost(9, COST:end) = [-30 -850 -10 -250  -5  -50 0 0];    % 10, 20, 30
%% put a load before gen in matrix
mpc.gen = [mpc.gen(8, :); mpc.gen(1:7, :); mpc.gen(9, :)];
mpc.gencost = [mpc.gencost(8, :); mpc.gencost(1:7, :); mpc.gencost(9, :)];
gg = find(~isload(mpc.gen));
ld = find( isload(mpc.gen));
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
orig_gc = mpc.gencost;

%%-----  single load zone, one scale factor  -----
dmd = 2;
t = 'all fixed loads (PQ) * 2 : ';
bus = scale_load(dmd, mpc.bus);
t_is(sum(bus(:, PD)), dmd*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), dmd*total.fixed.q, 8, [t 'total fixed Q']);
opt = struct('which', 'FIXED');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), dmd*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), dmd*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all fixed loads (P) * 2 : ';
opt = struct('pq', 'P');
bus = scale_load(dmd, mpc.bus, [], [], opt);
t_is(sum(bus(:, PD)), dmd*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
opt = struct('pq', 'P', 'which', 'FIXED');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), dmd*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all loads (PQ) * 2 : ';
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen);
t_is(sum(bus(:, PD)), dmd*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), dmd*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), dmd*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), dmd*total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all loads/costs (PQ) * 2 : ';
[bus, gen, gencost] = scale_load(dmd, mpc.bus, mpc.gen, [], [], mpc.gencost);
t_is(sum(bus(:, PD)), dmd*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), dmd*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), dmd*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), dmd*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end), 2*orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all loads/costs (PQ) * 2 : ';
opt = struct('cost', 1);
[bus, gen, gencost] = scale_load(dmd, mpc.bus, mpc.gen, [], opt, mpc.gencost);
t_is(sum(bus(:, PD)), dmd*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), dmd*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), dmd*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), dmd*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end), 2*orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all loads (PQ) * 2 : ';
opt = struct('cost', 0);
[bus, gen, gencost] = scale_load(dmd, mpc.bus, mpc.gen, [], opt, mpc.gencost);
t_is(sum(bus(:, PD)), dmd*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), dmd*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), dmd*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), dmd*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end),   orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all loads (P) * 2 : ';
opt = struct('pq', 'P');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), dmd*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all loads/costs (P) * 2 : ';
opt = struct('pq', 'P');
[bus, gen, gencost] = scale_load(dmd, mpc.bus, mpc.gen, [], opt, mpc.gencost);
t_is(sum(bus(:, PD)), dmd*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end), 2*orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all disp loads (PQ) * 2 : ';
opt = struct('which', 'DISPATCHABLE');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), dmd*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), dmd*total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all disp loads/costs (PQ) * 2 : ';
opt = struct('which', 'DISPATCHABLE');
[bus, gen, gencost] = scale_load(dmd, mpc.bus, mpc.gen, [], opt, mpc.gencost);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), dmd*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), dmd*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end), 2*orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all disp loads (P) * 2 : ';
opt = struct('pq', 'P', 'which', 'DISPATCHABLE');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all disp loads/costs (P) * 2 : ';
opt = struct('pq', 'P', 'which', 'DISPATCHABLE');
[bus, gen, gencost] = scale_load(dmd, mpc.bus, mpc.gen, [], opt, mpc.gencost);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end), 2*orig_gc(ld, COST:end), 8, [t 'gencost loads']);

%%-----  single load zone, one scale quantity  -----
dmd = 200;
t = 'all fixed loads (PQ) => total = 200 : ';
opt = struct('scale', 'QUANTITY');
bus = scale_load(dmd, mpc.bus, [], [], opt);
t_is(sum(bus(:, PD)), dmd, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), dmd/total.fixed.p*total.fixed.q, 8, [t 'total fixed Q']);
opt = struct('scale', 'QUANTITY', 'which', 'FIXED');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), dmd-total.disp.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), (dmd-total.disp.p)/total.fixed.p*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all fixed loads (P) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'pq', 'P');
bus = scale_load(dmd, mpc.bus, [], [], opt);
t_is(sum(bus(:, PD)), dmd, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
opt = struct('scale', 'QUANTITY', 'pq', 'P', 'which', 'FIXED');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), dmd-total.disp.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all loads (PQ) => total = 200 : ';
opt = struct('scale', 'QUANTITY');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), dmd/total.both.p*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), dmd/total.both.p*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd/total.both.p*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), dmd/total.both.p*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), dmd/total.both.p*total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all loads/costs (PQ) => total = 200 : ';
opt = struct('scale', 'QUANTITY');
[bus, gen, gencost] = scale_load(dmd, mpc.bus, mpc.gen, [], opt, mpc.gencost);
t_is(sum(bus(:, PD)), dmd/total.both.p*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), dmd/total.both.p*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd/total.both.p*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), dmd/total.both.p*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), dmd/total.both.p*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end), dmd/total.both.p*orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all loads (P) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'pq', 'P');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), dmd/total.both.p*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd/total.both.p*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all loads/costs (P) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'pq', 'P');
[bus, gen, gencost] = scale_load(dmd, mpc.bus, mpc.gen, [], opt, mpc.gencost);
t_is(sum(bus(:, PD)), dmd/total.both.p*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd/total.both.p*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end), dmd/total.both.p*orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all disp loads (PQ) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'which', 'DISPATCHABLE');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd-total.fixed.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), (dmd-total.fixed.p)/total.disp.p*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), (dmd-total.fixed.p)/total.disp.p*total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all disp loads/costs (PQ) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'which', 'DISPATCHABLE');
[bus, gen, gencost] = scale_load(dmd, mpc.bus, mpc.gen, [], opt, mpc.gencost);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd-total.fixed.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), (dmd-total.fixed.p)/total.disp.p*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), (dmd-total.fixed.p)/total.disp.p*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end), (dmd-total.fixed.p)/total.disp.p*orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all disp loads (P) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'pq', 'P', 'which', 'DISPATCHABLE');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd-total.fixed.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);

t = 'all disp loads/costs (P) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'pq', 'P', 'which', 'DISPATCHABLE');
[bus, gen, gencost] = scale_load(dmd, mpc.bus, mpc.gen, [], opt, mpc.gencost);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd-total.fixed.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end), (dmd-total.fixed.p)/total.disp.p*orig_gc(ld, COST:end), 8, [t 'gencost loads']);

%%-----  3 zones, area scale factors  -----
t = 'area fixed loads (PQ) * [3 2 1] : ';
dmd = [3 2 1];
bus = scale_load(dmd, mpc.bus);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), dmd(k)*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
end
opt = struct('which', 'FIXED');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), dmd(k)*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

t = 'area fixed loads (P) * [3 2 1] : ';
dmd = [3 2 1];
opt = struct('pq', 'P');
bus = scale_load(dmd, mpc.bus, [], [], opt);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
end
opt = struct('pq', 'P', 'which', 'FIXED');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

t = 'all area loads (PQ) * [3 2 1] : ';
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), dmd(k)*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), dmd(k)*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), dmd(k)*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

t = 'all area loads/costs (PQ) * [3 2 1] : ';
[bus, gen, gencost] = scale_load(dmd, mpc.bus, mpc.gen, [], [], mpc.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), dmd(k)*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), dmd(k)*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), dmd(k)*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    t_is(gencost(lda{k}, COST:end), dmd(k)*orig_gc(lda{k}, COST:end), 8, sprintf('%s area %d gencost loads', t, k));
end
t_is(gencost(gg, COST:end), orig_gc(gg, COST:end), 8, sprintf('%s gencost gens', t));

t = 'all area loads (P) * [3 2 1] : ';
opt = struct('pq', 'P');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

t = 'all area loads/costs (P) * [3 2 1] : ';
opt = struct('pq', 'P');
[bus, gen, gencost] = scale_load(dmd, mpc.bus, mpc.gen, [], opt, mpc.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    t_is(gencost(lda{k}, COST:end), dmd(k)*orig_gc(lda{k}, COST:end), 8, sprintf('%s area %d gencost loads', t, k));
end
t_is(gencost(gg, COST:end), orig_gc(gg, COST:end), 8, sprintf('%s gencost gens', t));

t = 'area disp loads (PQ) * [3 2 1] : ';
opt = struct('which', 'DISPATCHABLE');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), dmd(k)*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), dmd(k)*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

t = 'area disp loads/costs (PQ) * [3 2 1] : ';
opt = struct('which', 'DISPATCHABLE');
[bus, gen, gencost] = scale_load(dmd, mpc.bus, mpc.gen, [], opt, mpc.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), dmd(k)*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), dmd(k)*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    t_is(gencost(lda{k}, COST:end), dmd(k)*orig_gc(lda{k}, COST:end), 8, sprintf('%s area %d gencost loads', t, k));
end
t_is(gencost(gg, COST:end), orig_gc(gg, COST:end), 8, sprintf('%s gencost gens', t));

t = 'area disp loads (P) * [3 2 1] : ';
opt = struct('pq', 'P', 'which', 'DISPATCHABLE');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

t = 'area disp loads/costs (P) * [3 2 1] : ';
opt = struct('pq', 'P', 'which', 'DISPATCHABLE');
[bus, gen, gencost] = scale_load(dmd, mpc.bus, mpc.gen, [], opt, mpc.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    t_is(gencost(lda{k}, COST:end), dmd(k)*orig_gc(lda{k}, COST:end), 8, sprintf('%s area %d gencost loads', t, k));
end
t_is(gencost(gg, COST:end), orig_gc(gg, COST:end), 8, sprintf('%s gencost gens', t));

%%-----  3 zones, area scale quantities  -----
t = 'area fixed loads (PQ) => total = [100 80 60] : ';
dmd = [100 80 60];
opt = struct('scale', 'QUANTITY');
bus = scale_load(dmd, mpc.bus, [], [], opt);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k), 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), dmd(k)/area(k).fixed.p*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
end
opt = struct('scale', 'QUANTITY', 'which', 'FIXED');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)-area(k).disp.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), (dmd(k)-area(k).disp.p)/area(k).fixed.p*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

t = 'area fixed loads (P) => total = [100 80 60] : ';
dmd = [100 80 60];
opt = struct('scale', 'QUANTITY', 'pq', 'P');
bus = scale_load(dmd, mpc.bus, [], [], opt);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k), 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
end
opt = struct('scale', 'QUANTITY', 'pq', 'P', 'which', 'FIXED');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)-area(k).disp.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

t = 'all area loads (PQ) => total = [100 80 60] : ';
opt = struct('scale', 'QUANTITY');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)/area(k).both.p*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), dmd(k)/area(k).both.p*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)/area(k).both.p*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), dmd(k)/area(k).both.p*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), dmd(k)/area(k).both.p*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

t = 'all area loads (P) => total = [100 80 60] : ';
opt = struct('scale', 'QUANTITY', 'pq', 'P');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)/area(k).both.p*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)/area(k).both.p*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

t = 'area disp loads (PQ) => total = [100 80 60] : throws expected exception';
dmd = [100 80 60];
opt = struct('scale', 'QUANTITY', 'which', 'DISPATCHABLE');
err = 0;
try
    [bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
catch
    [msg, id] = lasterr;
    expected = 'scale_load: impossible to make zone 2 load equal 80 by scaling non-existent dispatchable load';
    if ~isempty(findstr(expected, msg))
        err = 1;
    end
end
t_ok(err, t);

t = 'area disp loads (PQ) => total = [100 74.3941 60] : ';
dmd = [100 area(2).fixed.p 60];
opt = struct('scale', 'QUANTITY', 'which', 'DISPATCHABLE');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)-area(k).fixed.p, 8, sprintf('%s area %d disp P', t, k));
    if k == 2
        t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    else
        t_is(-sum(gen(lda{k}, QMIN)), (dmd(k)-area(k).fixed.p)/area(k).disp.p*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(gen(lda{k}, QMAX)), (dmd(k)-area(k).fixed.p)/area(k).disp.p*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

t = 'area disp loads (P) => total = [100 74.3941 60] : ';
opt = struct('scale', 'QUANTITY', 'pq', 'P', 'which', 'DISPATCHABLE');
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, [], opt);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)-area(k).fixed.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

%%-----  explict single load zone  -----
t = 'explicit single load zone';
load_zone = zeros(1, size(mpc.bus, 1));
load_zone([3 4]) = 1;
dmd = 2;
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, load_zone);
Pd = mpc.bus(:, PD);
Pd([3 4]) = dmd * Pd([3 4]);
t_is( bus(:, PD), Pd, 8, t);

%%-----  explict multiple load zone  -----
t = 'explicit multiple load zone';
load_zone = zeros(1, size(mpc.bus, 1));
load_zone([3 4]) = 1;
load_zone([7 8]) = 2;
dmd = [2 0.5];
[bus, gen] = scale_load(dmd, mpc.bus, mpc.gen, load_zone);
Pd = mpc.bus(:, PD);
Pd([3 4]) = dmd(1) * Pd([3 4]);
Pd([7 8]) = dmd(2) * Pd([7 8]);
t_is( bus(:, PD), Pd, 8, t);



%%------------------------------------------------------------------
%% mostly same tests below, but with MPC input/output

%%-----  single load zone, one scale factor  -----
dmd = 2;
t = 'all loads (PQ) * 2 : ';
opt = struct('cost', 0);
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
t_is(sum(bus(:, PD)), dmd*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), dmd*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), dmd*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), dmd*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end),   orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all loads/costs (PQ) * 2 : ';
mpc1 = scale_load(dmd, mpc);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
t_is(sum(bus(:, PD)), dmd*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), dmd*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), dmd*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), dmd*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end), 2*orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all loads (P) * 2 : ';
opt = struct('pq', 'P', 'cost', 0);
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
t_is(sum(bus(:, PD)), dmd*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end),   orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all loads/costs (P) * 2 : ';
opt = struct('pq', 'P');
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
t_is(sum(bus(:, PD)), dmd*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end), 2*orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all disp loads (PQ) * 2 : ';
opt = struct('which', 'DISPATCHABLE', 'cost', 0);
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), dmd*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), dmd*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end),   orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all disp loads/costs (PQ) * 2 : ';
opt = struct('which', 'DISPATCHABLE');
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), dmd*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), dmd*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end), 2*orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all disp loads (P) * 2 : ';
opt = struct('pq', 'P', 'which', 'DISPATCHABLE', 'cost', 0);
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end),   orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all disp loads/costs (P) * 2 : ';
opt = struct('pq', 'P', 'which', 'DISPATCHABLE', 'cost', 1);
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end), 2*orig_gc(ld, COST:end), 8, [t 'gencost loads']);

%%-----  single load zone, one scale quantity  -----
dmd = 200;
t = 'all fixed loads (PQ) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'which', 'FIXED');
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
t_is(sum(bus(:, PD)), dmd-total.disp.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), (dmd-total.disp.p)/total.fixed.p*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end),   orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all fixed loads (P) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'pq', 'P', 'which', 'FIXED', 'cost', 1);
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
t_is(sum(bus(:, PD)), dmd-total.disp.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end),   orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all loads (PQ) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'cost', 0);
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
t_is(sum(bus(:, PD)), dmd/total.both.p*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), dmd/total.both.p*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd/total.both.p*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), dmd/total.both.p*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), dmd/total.both.p*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end),   orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all loads/costs (PQ) => total = 200 : ';
opt = struct('scale', 'QUANTITY');
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
t_is(sum(bus(:, PD)), dmd/total.both.p*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), dmd/total.both.p*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd/total.both.p*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), dmd/total.both.p*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), dmd/total.both.p*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end), dmd/total.both.p*orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all loads (P) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'pq', 'P', 'cost', 0);
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
t_is(sum(bus(:, PD)), dmd/total.both.p*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd/total.both.p*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end),   orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all loads/costs (P) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'pq', 'P');
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
t_is(sum(bus(:, PD)), dmd/total.both.p*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd/total.both.p*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end), dmd/total.both.p*orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all disp loads (PQ) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'which', 'DISPATCHABLE', 'cost', 0);
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd-total.fixed.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), (dmd-total.fixed.p)/total.disp.p*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), (dmd-total.fixed.p)/total.disp.p*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end),   orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all disp loads/costs (PQ) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'which', 'DISPATCHABLE');
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd-total.fixed.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), (dmd-total.fixed.p)/total.disp.p*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), (dmd-total.fixed.p)/total.disp.p*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end), (dmd-total.fixed.p)/total.disp.p*orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all disp loads (P) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'pq', 'P', 'which', 'DISPATCHABLE', 'cost', 0);
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd-total.fixed.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end),   orig_gc(ld, COST:end), 8, [t 'gencost loads']);

t = 'all disp loads/costs (P) => total = 200 : ';
opt = struct('scale', 'QUANTITY', 'pq', 'P', 'which', 'DISPATCHABLE');
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
t_is(sum(bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(gen(ld, PMIN)), dmd-total.fixed.p, 8, [t 'total disp P']);
t_is(-sum(gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(gencost(gg, COST:end),   orig_gc(gg, COST:end), 8, [t 'gencost gens']);
t_is(gencost(ld, COST:end), (dmd-total.fixed.p)/total.disp.p*orig_gc(ld, COST:end), 8, [t 'gencost loads']);

%%-----  3 zones, area scale factors  -----
t = 'area fixed loads (PQ) * [3 2 1] : ';
dmd = [3 2 1];
opt = struct('which', 'FIXED', 'cost', 1);
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), dmd(k)*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    t_is(gencost(lda{k}, COST:end), orig_gc(lda{k}, COST:end), 8, sprintf('%s area %d gencost loads', t, k));
end
t_is(gencost(gg, COST:end), orig_gc(gg, COST:end), 8, sprintf('%s gencost gens', t));

t = 'area fixed loads (P) * [3 2 1] : ';
dmd = [3 2 1];
opt = struct('pq', 'P', 'which', 'FIXED', 'cost', 1);
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    t_is(gencost(lda{k}, COST:end), orig_gc(lda{k}, COST:end), 8, sprintf('%s area %d gencost loads', t, k));
end
t_is(gencost(gg, COST:end), orig_gc(gg, COST:end), 8, sprintf('%s gencost gens', t));

t = 'all area loads (PQ) * [3 2 1] : ';
opt = struct('cost', 0);
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), dmd(k)*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), dmd(k)*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), dmd(k)*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    t_is(gencost(lda{k}, COST:end), orig_gc(lda{k}, COST:end), 8, sprintf('%s area %d gencost loads', t, k));
end
t_is(gencost(gg, COST:end), orig_gc(gg, COST:end), 8, sprintf('%s gencost gens', t));

t = 'all area loads/costs (PQ) * [3 2 1] : ';
mpc1 = scale_load(dmd, mpc);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), dmd(k)*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), dmd(k)*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), dmd(k)*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    t_is(gencost(lda{k}, COST:end), dmd(k)*orig_gc(lda{k}, COST:end), 8, sprintf('%s area %d gencost loads', t, k));
end
t_is(gencost(gg, COST:end), orig_gc(gg, COST:end), 8, sprintf('%s gencost gens', t));

t = 'all area loads (P) * [3 2 1] : ';
opt = struct('pq', 'P', 'cost', 0);
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    t_is(gencost(lda{k}, COST:end), orig_gc(lda{k}, COST:end), 8, sprintf('%s area %d gencost loads', t, k));
end
t_is(gencost(gg, COST:end), orig_gc(gg, COST:end), 8, sprintf('%s gencost gens', t));

t = 'all area loads/costs (P) * [3 2 1] : ';
opt = struct('pq', 'P');
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    t_is(gencost(lda{k}, COST:end), dmd(k)*orig_gc(lda{k}, COST:end), 8, sprintf('%s area %d gencost loads', t, k));
end
t_is(gencost(gg, COST:end), orig_gc(gg, COST:end), 8, sprintf('%s gencost gens', t));

t = 'area disp loads (PQ) * [3 2 1] : ';
opt = struct('which', 'DISPATCHABLE', 'cost', 0);
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), dmd(k)*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), dmd(k)*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    t_is(gencost(lda{k}, COST:end), orig_gc(lda{k}, COST:end), 8, sprintf('%s area %d gencost loads', t, k));
end
t_is(gencost(gg, COST:end), orig_gc(gg, COST:end), 8, sprintf('%s gencost gens', t));

t = 'area disp loads/costs (PQ) * [3 2 1] : ';
opt = struct('which', 'DISPATCHABLE');
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), dmd(k)*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), dmd(k)*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    t_is(gencost(lda{k}, COST:end), dmd(k)*orig_gc(lda{k}, COST:end), 8, sprintf('%s area %d gencost loads', t, k));
end
t_is(gencost(gg, COST:end), orig_gc(gg, COST:end), 8, sprintf('%s gencost gens', t));

t = 'area disp loads (P) * [3 2 1] : ';
opt = struct('pq', 'P', 'which', 'DISPATCHABLE', 'cost', 0);
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    t_is(gencost(lda{k}, COST:end), orig_gc(lda{k}, COST:end), 8, sprintf('%s area %d gencost loads', t, k));
end
t_is(gencost(gg, COST:end), orig_gc(gg, COST:end), 8, sprintf('%s gencost gens', t));

t = 'area disp loads/costs (P) * [3 2 1] : ';
opt = struct('pq', 'P', 'which', 'DISPATCHABLE');
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    t_is(gencost(lda{k}, COST:end), dmd(k)*orig_gc(lda{k}, COST:end), 8, sprintf('%s area %d gencost loads', t, k));
end
t_is(gencost(gg, COST:end), orig_gc(gg, COST:end), 8, sprintf('%s gencost gens', t));

%%-----  3 zones, area scale quantities  -----
t = 'area fixed loads (PQ) => total = [100 80 60] : ';
dmd = [100 80 60];
opt = struct('scale', 'QUANTITY', 'which', 'FIXED');
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)-area(k).disp.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), (dmd(k)-area(k).disp.p)/area(k).fixed.p*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    t_is(gencost(lda{k}, COST:end), orig_gc(lda{k}, COST:end), 8, sprintf('%s area %d gencost loads', t, k));
end
t_is(gencost(gg, COST:end), orig_gc(gg, COST:end), 8, sprintf('%s gencost gens', t));

t = 'area fixed loads (P) => total = [100 80 60] : ';
dmd = [100 80 60];
opt = struct('scale', 'QUANTITY', 'pq', 'P', 'which', 'FIXED', 'cost', 1);
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)-area(k).disp.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    t_is(gencost(lda{k}, COST:end), orig_gc(lda{k}, COST:end), 8, sprintf('%s area %d gencost loads', t, k));
end
t_is(gencost(gg, COST:end), orig_gc(gg, COST:end), 8, sprintf('%s gencost gens', t));

t = 'all area loads (PQ) => total = [100 80 60] : ';
opt = struct('scale', 'QUANTITY', 'cost', 0);
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)/area(k).both.p*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), dmd(k)/area(k).both.p*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)/area(k).both.p*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), dmd(k)/area(k).both.p*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), dmd(k)/area(k).both.p*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    t_is(gencost(lda{k}, COST:end), orig_gc(lda{k}, COST:end), 8, sprintf('%s area %d gencost loads', t, k));
end
t_is(gencost(gg, COST:end), orig_gc(gg, COST:end), 8, sprintf('%s gencost gens', t));

t = 'all area loads/costs (P) => total = [100 80 60] : ';
opt = struct('scale', 'QUANTITY', 'pq', 'P');
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), dmd(k)/area(k).both.p*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)/area(k).both.p*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    t_is(gencost(lda{k}, COST:end), dmd(k)/area(k).both.p*orig_gc(lda{k}, COST:end), 8, sprintf('%s area %d gencost loads', t, k));
end
t_is(gencost(gg, COST:end), orig_gc(gg, COST:end), 8, sprintf('%s gencost gens', t));

t = 'area disp loads (PQ) => total = [100 80 60] : throws expected exception';
dmd = [100 80 60];
opt = struct('scale', 'QUANTITY', 'which', 'DISPATCHABLE');
err = 0;
try
    mpc1 = scale_load(dmd, mpc, [], opt);
catch
    [msg, id] = lasterr;
    expected = 'scale_load: impossible to make zone 2 load equal 80 by scaling non-existent dispatchable load';
    if ~isempty(findstr(expected, msg))
        err = 1;
    end
end
t_ok(err, t);

t = 'area disp loads (PQ) => total = [100 74.3941 60] : ';
dmd = [100 area(2).fixed.p 60];
opt = struct('scale', 'QUANTITY', 'which', 'DISPATCHABLE');
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)-area(k).fixed.p, 8, sprintf('%s area %d disp P', t, k));
    if k == 2
        t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    else
        t_is(-sum(gen(lda{k}, QMIN)), (dmd(k)-area(k).fixed.p)/area(k).disp.p*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(gen(lda{k}, QMAX)), (dmd(k)-area(k).fixed.p)/area(k).disp.p*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

t = 'area disp loads (P) => total = [100 74.3941 60] : ';
opt = struct('scale', 'QUANTITY', 'pq', 'P', 'which', 'DISPATCHABLE');
mpc1 = scale_load(dmd, mpc, [], opt);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
for k = 1:length(dmd)
    t_is(sum(bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
    t_is(sum(bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
    t_is(-sum(gen(lda{k}, PMIN)), dmd(k)-area(k).fixed.p, 8, sprintf('%s area %d disp P', t, k));
    t_is(-sum(gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
    t_is(-sum(gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
end

%%-----  explict single load zone  -----
t = 'explicit single load zone';
load_zone = zeros(1, size(mpc.bus, 1));
load_zone([3 4]) = 1;
dmd = 2;
mpc1 = scale_load(dmd, mpc, load_zone);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
Pd = mpc.bus(:, PD);
Pd([3 4]) = dmd * Pd([3 4]);
t_is( bus(:, PD), Pd, 8, t);

%%-----  explict multiple load zone  -----
t = 'explicit multiple load zone';
load_zone = zeros(1, size(mpc.bus, 1));
load_zone([3 4]) = 1;
load_zone([7 8]) = 2;
dmd = [2 0.5];
mpc1 = scale_load(dmd, mpc, load_zone);
[bus, gen, gencost] = deal(mpc1.bus, mpc1.gen, mpc1.gencost);
Pd = mpc.bus(:, PD);
Pd([3 4]) = dmd(1) * Pd([3 4]);
Pd([7 8]) = dmd(2) * Pd([7 8]);
t_is( bus(:, PD), Pd, 8, t);

t_end;
