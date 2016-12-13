function t_total_load(quiet)
%T_TOTAL_LOAD  Tests for code in TOTAL_LOAD.

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

n_tests = 250;

t_begin(n_tests, quiet);

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

mpc = loadcase('t_auction_case');
mpc.gen(8, GEN_BUS) = 2;    %% multiple d. loads per area, same bus as gen
mpc.gen(8, [PG QG QMIN QMAX]) = [-20 2 0 3 ];
%% put load before gen in matrix
mpc.gen = [mpc.gen(8, :); mpc.gen(1:7, :); mpc.gen(9, :)];
%% add an isolated bus (with both kinds of load)
k = 12;
mpc.bus = [mpc.bus(1:k, :); mpc.bus(k, :); mpc.bus(k+1:end, :)];
mpc.bus(k, BUS_I) = mpc.bus(k, BUS_I)*10;
mpc.bus(k, BUS_TYPE) = NONE;
mpc.gen = [mpc.gen; mpc.gen(end, :)];
mpc.gen(end, GEN_BUS) = mpc.bus(k, BUS_I);

ld = find(isload(mpc.gen(1:end-1, :)));
for k = 1:3
    a{k} = find(mpc.bus(:, BUS_AREA) == k & mpc.bus(:, BUS_TYPE) ~= NONE); %% buses in area k
    [junk, tmp, junk2] = intersect(mpc.gen(ld, GEN_BUS), a{k});
    lda{k} = ld(tmp);                       %% disp loads in area k
end
for k = 1:3
    area(k).fixed.p = sum(mpc.bus(a{k}, PD));
    area(k).fixed.q = sum(mpc.bus(a{k}, QD));
    area(k).disp.pnom = -sum(mpc.gen(lda{k}, PMIN));
    area(k).disp.p = -sum(mpc.gen(lda{k}, PG));
    area(k).disp.qmin = -sum(mpc.gen(lda{k}, QMIN));
    area(k).disp.qmax = -sum(mpc.gen(lda{k}, QMAX));
    area(k).disp.qnom = area(k).disp.qmin + area(k).disp.qmax;
    area(k).disp.q = -sum(mpc.gen(lda{k}, QG));
    area(k).both.pnom = area(k).fixed.p + area(k).disp.pnom;
    area(k).both.qnom = area(k).fixed.q + area(k).disp.qnom;
    area(k).both.p = area(k).fixed.p + area(k).disp.p;
    area(k).both.q = area(k).fixed.q + area(k).disp.q;
end
k = find(mpc.bus(:, BUS_TYPE) ~= NONE);
total.fixed.p = sum(mpc.bus(k, PD));
total.fixed.q = sum(mpc.bus(k, QD));
total.disp.p = -sum(mpc.gen(ld, PG));
total.disp.pnom = -sum(mpc.gen(ld, PMIN));
total.disp.qmin = -sum(mpc.gen(ld, QMIN));
total.disp.qmax = -sum(mpc.gen(ld, QMAX));
total.disp.q = -sum(mpc.gen(ld, QG));
total.disp.qnom = total.disp.qmin + total.disp.qmax;
total.both.pnom = total.fixed.p + total.disp.pnom;
total.both.qnom = total.fixed.q + total.disp.qnom;
total.both.p = total.fixed.p + total.disp.p;
total.both.q = total.fixed.q + total.disp.q;

%%-----  all load  -----
t = '      Pd = total_load(bus) : ';
Pd = total_load(mpc.bus);
t_is(Pd, [area(1).fixed.p; area(2).fixed.p; area(3).fixed.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(bus) : ';
[Pd, Qd] = total_load(mpc.bus);
t_is(Pd, [area(1).fixed.p; area(2).fixed.p; area(3).fixed.p], 12, [t 'Pd']);
t_is(Qd, [area(1).fixed.q; area(2).fixed.q; area(3).fixed.q], 12, [t 'Qd']);

t = '      Pd = total_load(bus, gen) : ';
Pd = total_load(mpc.bus, mpc.gen);
t_is(Pd, [area(1).both.p; area(2).both.p; area(3).both.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(bus, gen) : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen);
t_is(Pd, [area(1).both.p; area(2).both.p; area(3).both.p], 12, [t 'Pd']);
t_is(Qd, [area(1).both.q; area(2).both.q; area(3).both.q], 12, [t 'Qd']);

t = '      Pd = total_load(bus, [], ''all'') : ';
Pd = total_load(mpc.bus, [], 'all');
t_is(Pd, total.fixed.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(bus, [], ''all'') : ';
[Pd, Qd] = total_load(mpc.bus, [], 'all');
t_is(Pd, total.fixed.p, 12, [t 'Pd']);
t_is(Qd, total.fixed.q, 12, [t 'Qd']);

t = '      Pd = total_load(b, g, ''all'') : ';
Pd = total_load(mpc.bus, mpc.gen, 'all');
t_is(Pd, total.both.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, ''all'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, 'all');
t_is(Pd, total.both.p, 12, [t 'Pd']);
t_is(Qd, total.both.q, 12, [t 'Qd']);

t = '      Pd = total_load(b, g, ''all'', .type = ''BOTH'') : ';
opt = struct('type', 'BOTH');
Pd = total_load(mpc.bus, mpc.gen, 'all', opt);
t_is(Pd, total.both.p, 12, [t 'Pd']);

t = '      Pd = total_load(b, g, ''all'', .type = ''BOTH'' .nominal = 1) : ';
opt = struct('type', 'BOTH', 'nominal', 1);
Pd = total_load(mpc.bus, mpc.gen, 'all', opt);
t_is(Pd, total.both.pnom, 12, [t 'Pd']);

t = '      Pd = total_load(b, g, ''all'', .type = ''BOTH'' .nominal = 0) : ';
opt = struct('type', 'BOTH', 'nominal', 0);
Pd = total_load(mpc.bus, mpc.gen, 'all', opt);
t_is(Pd, total.both.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, ''all'', .type = ''BOTH'') : ';
opt = struct('type', 'BOTH');
[Pd, Qd] = total_load(mpc.bus, mpc.gen, 'all', opt);
t_is(Pd, total.both.p, 12, [t 'Pd']);
t_is(Qd, total.both.q, 12, [t 'Qd']);

t = '      Pd = total_load(b, g, ''all'', .type = ''FIXED'') : ';
opt = struct('type', 'FIXED');
Pd = total_load(mpc.bus, mpc.gen, 'all', opt);
t_is(Pd, total.fixed.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, ''all'', .type = ''FIXED'') : ';
opt = struct('type', 'FIXED');
[Pd, Qd] = total_load(mpc.bus, mpc.gen, 'all', opt);
t_is(Pd, total.fixed.p, 12, [t 'Pd']);
t_is(Qd, total.fixed.q, 12, [t 'Qd']);

t = '      Pd = total_load(b, g, ''all'', .type = ''DISPATCHABLE'') : ';
opt = struct('type', 'DISPATCHABLE');
Pd = total_load(mpc.bus, mpc.gen, 'all', opt);
t_is(Pd, total.disp.p, 12, [t 'Pd']);

t = '      Pd = total_load(b, g, ''all'', .type = ''DISP'' .nominal = 1) : ';
opt = struct('type', 'DISP', 'nominal', 1);
Pd = total_load(mpc.bus, mpc.gen, 'all', opt);
t_is(Pd, total.disp.pnom, 12, [t 'Pd']);

t = '      Pd = total_load(b, g, ''all'', .type = ''DISP'' .nominal = 0) : ';
opt = struct('type', 'DISP', 'nominal', 0);
Pd = total_load(mpc.bus, mpc.gen, 'all', opt);
t_is(Pd, total.disp.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, ''all'', .type = ''DISPATCHABLE'') : ';
opt = struct('type', 'DISPATCHABLE');
[Pd, Qd] = total_load(mpc.bus, mpc.gen, 'all', opt);
t_is(Pd, total.disp.p, 12, [t 'Pd']);
t_is(Qd, total.disp.q, 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(b, g, ''all'', .type = ''DISP'' .nominal = 1) : ';
opt = struct('type', 'DISP', 'nominal', 1);
[Pd, Qd] = total_load(mpc.bus, mpc.gen, 'all', opt);
t_is(Pd, total.disp.pnom, 12, [t 'Pd']);
t_is(Qd, total.disp.qnom, 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(b, g, ''all'', .type = ''DISP'' .nominal = 0) : ';
opt = struct('type', 'DISP', 'nominal', 0);
[Pd, Qd] = total_load(mpc.bus, mpc.gen, 'all', opt);
t_is(Pd, total.disp.p, 12, [t 'Pd']);
t_is(Qd, total.disp.q, 12, [t 'Qd']);

t = '      Pd = total_load(b, g, [], .type = ''BOTH'') : ';
opt = struct('type', 'BOTH');
Pd = total_load(mpc.bus, mpc.gen, [], opt);
t_is(Pd, [area(1).both.p; area(2).both.p; area(3).both.p], 12, [t 'Pd']);

t = '      Pd = total_load(b, g, [], .type = ''BOTH'' .nominal = 1) : ';
opt = struct('type', 'BOTH', 'nominal', 1);
Pd = total_load(mpc.bus, mpc.gen, [], opt);
t_is(Pd, [area(1).both.pnom; area(2).both.pnom; area(3).both.pnom], 12, [t 'Pd']);

t = '      Pd = total_load(b, g, [], .type = ''BOTH'' .nominal = 0) : ';
opt = struct('type', 'BOTH', 'nominal', 0);
Pd = total_load(mpc.bus, mpc.gen, [], opt);
t_is(Pd, [area(1).both.p; area(2).both.p; area(3).both.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, [], .type = ''BOTH'') : ';
opt = struct('type', 'BOTH');
[Pd, Qd] = total_load(mpc.bus, mpc.gen, [], opt);
t_is(Pd, [area(1).both.p; area(2).both.p; area(3).both.p], 12, [t 'Pd']);
t_is(Qd, [area(1).both.q; area(2).both.q; area(3).both.q], 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(b, g, [], .type = ''BOTH'' .nominal = 1) : ';
opt = struct('type', 'BOTH', 'nominal', 1);
[Pd, Qd] = total_load(mpc.bus, mpc.gen, [], opt);
t_is(Pd, [area(1).both.pnom; area(2).both.pnom; area(3).both.pnom], 12, [t 'Pd']);
t_is(Qd, [area(1).both.qnom; area(2).both.qnom; area(3).both.qnom], 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(b, g, [], .type = ''BOTH'' .nominal = 0) : ';
opt = struct('type', 'BOTH', 'nominal', 0);
[Pd, Qd] = total_load(mpc.bus, mpc.gen, [], opt);
t_is(Pd, [area(1).both.p; area(2).both.p; area(3).both.p], 12, [t 'Pd']);
t_is(Qd, [area(1).both.q; area(2).both.q; area(3).both.q], 12, [t 'Qd']);

t = '      Pd = total_load(b, g, [], .type = ''FIXED'') : ';
opt = struct('type', 'FIXED');
Pd = total_load(mpc.bus, mpc.gen, [], opt);
t_is(Pd, [area(1).fixed.p; area(2).fixed.p; area(3).fixed.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, [], .type = ''FIXED'') : ';
opt = struct('type', 'FIXED');
[Pd, Qd] = total_load(mpc.bus, mpc.gen, [], opt);
t_is(Pd, [area(1).fixed.p; area(2).fixed.p; area(3).fixed.p], 12, [t 'Pd']);
t_is(Qd, [area(1).fixed.q; area(2).fixed.q; area(3).fixed.q], 12, [t 'Qd']);

t = '      Pd = total_load(b, g, [], .type = ''DISPATCHABLE'') : ';
opt = struct('type', 'DISPATCHABLE');
Pd = total_load(mpc.bus, mpc.gen, [], opt);
t_is(Pd, [area(1).disp.p; area(2).disp.p; area(3).disp.p], 12, [t 'Pd']);

t = '      Pd = total_load(b, g, [], .type = ''DISP'' .nominal = 1) : ';
opt = struct('type', 'DISP', 'nominal', 1);
Pd = total_load(mpc.bus, mpc.gen, [], opt);
t_is(Pd, [area(1).disp.pnom; area(2).disp.pnom; area(3).disp.pnom], 12, [t 'Pd']);

t = '      Pd = total_load(b, g, [], .type = ''DISP'' .nominal = 0) : ';
opt = struct('type', 'DISP', 'nominal', 0);
Pd = total_load(mpc.bus, mpc.gen, [], opt);
t_is(Pd, [area(1).disp.p; area(2).disp.p; area(3).disp.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, [], .type = ''DISPATCHABLE'') : ';
opt = struct('type', 'DISPATCHABLE');
[Pd, Qd] = total_load(mpc.bus, mpc.gen, [], opt);
t_is(Pd, [area(1).disp.p; area(2).disp.p; area(3).disp.p], 12, [t 'Pd']);
t_is(Qd, [area(1).disp.q; area(2).disp.q; area(3).disp.q], 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(b, g, [], .type = ''DISP'' .nominal = 1) : ';
opt = struct('type', 'DISP', 'nominal', 1);
[Pd, Qd] = total_load(mpc.bus, mpc.gen, [], opt);
t_is(Pd, [area(1).disp.pnom; area(2).disp.pnom; area(3).disp.pnom], 12, [t 'Pd']);
t_is(Qd, [area(1).disp.qnom; area(2).disp.qnom; area(3).disp.qnom], 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(b, g, [], .type = ''DISP'' .nominal = 0) : ';
opt = struct('type', 'DISP', 'nominal', 0);
[Pd, Qd] = total_load(mpc.bus, mpc.gen, [], opt);
t_is(Pd, [area(1).disp.p; area(2).disp.p; area(3).disp.p], 12, [t 'Pd']);
t_is(Qd, [area(1).disp.q; area(2).disp.q; area(3).disp.q], 12, [t 'Qd']);

%%-----  explicit single load zone  -----
nb = size(mpc.bus, 1);
load_zone = zeros(nb, 1);
k = find(mpc.bus(:, BUS_AREA) == 1);    %% area 1
load_zone(k) = 1;
t = '      Pd = total_load(b, g, ld_zone1, .type = ''BOTH'') : ';
opt = struct('type', 'BOTH');
Pd = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, area(1).both.p, 12, [t 'Pd']);

t = '      Pd = total_load(b, g, ld_zone1, .type = ''BOTH'' .nominal = 1) : ';
opt = struct('type', 'BOTH', 'nominal', 1);
Pd = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, area(1).both.pnom, 12, [t 'Pd']);

t = '      Pd = total_load(b, g, ld_zone1, .type = ''BOTH'' .nominal = 0) : ';
opt = struct('type', 'BOTH', 'nominal', 0);
Pd = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, area(1).both.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone1, .type = ''BOTH'') : ';
opt = struct('type', 'BOTH');
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, area(1).both.p, 12, [t 'Pd']);
t_is(Qd, area(1).both.q, 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone1, .type = ''BOTH'' .nominal = 1) : ';
opt = struct('type', 'BOTH', 'nominal', 1);
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, area(1).both.pnom, 12, [t 'Pd']);
t_is(Qd, area(1).both.qnom, 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone1, .type = ''BOTH'' .nominal = 0) : ';
opt = struct('type', 'BOTH', 'nominal', 0);
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, area(1).both.p, 12, [t 'Pd']);
t_is(Qd, area(1).both.q, 12, [t 'Qd']);

t = '      Pd = total_load(b, g, ld_zone1, .type = ''FIXED'') : ';
opt = struct('type', 'FIXED');
Pd = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, area(1).fixed.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone1, .type = ''FIXED'') : ';
opt = struct('type', 'FIXED');
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, area(1).fixed.p, 12, [t 'Pd']);
t_is(Qd, area(1).fixed.q, 12, [t 'Qd']);

t = '      Pd = total_load(b, g, ld_zone1, .type = ''DISPATCHABLE'') : ';
opt = struct('type', 'DISPATCHABLE');
Pd = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, area(1).disp.p, 12, [t 'Pd']);

t = '      Pd = total_load(b, g, ld_zone1, .type = ''DISP'' .nominal = 1) : ';
opt = struct('type', 'DISP', 'nominal', 1);
Pd = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, area(1).disp.pnom, 12, [t 'Pd']);

t = '      Pd = total_load(b, g, ld_zone1, .type = ''DISP'' .nominal = 0) : ';
opt = struct('type', 'DISP', 'nominal', 0);
Pd = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, area(1).disp.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone1, .type = ''DISPATCHABLE'') : ';
opt = struct('type', 'DISPATCHABLE');
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, area(1).disp.p, 12, [t 'Pd']);
t_is(Qd, area(1).disp.q, 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone1, .type = ''DISP'' .nominal = 1) : ';
opt = struct('type', 'DISP', 'nominal', 1);
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, area(1).disp.pnom, 12, [t 'Pd']);
t_is(Qd, area(1).disp.qnom, 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone1, .type = ''DISP'' .nominal = 0) : ';
opt = struct('type', 'DISP', 'nominal', 0);
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, area(1).disp.p, 12, [t 'Pd']);
t_is(Qd, area(1).disp.q, 12, [t 'Qd']);

%%-----  explicit multiple load zone  -----
load_zone = zeros(nb, 1);
k = find(mpc.bus(:, BUS_AREA) == 3);    %% area 3
load_zone(k) = 1;
k = find(mpc.bus(:, BUS_AREA) == 1);    %% area 1
load_zone(k) = 2;
t = '      Pd = total_load(b, g, ld_zone2, .type = ''BOTH'') : ';
opt = struct('type', 'BOTH');
Pd = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, [area(3).both.p; area(1).both.p], 12, [t 'Pd']);

t = '      Pd = total_load(b, g, ld_zone2, .type = ''BOTH'' .nominal = 1) : ';
opt = struct('type', 'BOTH', 'nominal', 1);
Pd = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, [area(3).both.pnom; area(1).both.pnom], 12, [t 'Pd']);

t = '      Pd = total_load(b, g, ld_zone2, .type = ''BOTH'' .nominal = 0) : ';
opt = struct('type', 'BOTH', 'nominal', 0);
Pd = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, [area(3).both.p; area(1).both.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone2, .type = ''BOTH'') : ';
opt = struct('type', 'BOTH');
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, [area(3).both.p; area(1).both.p], 12, [t 'Pd']);
t_is(Qd, [area(3).both.q; area(1).both.q], 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone2, .type = ''BOTH'' .nominal = 1) : ';
opt = struct('type', 'BOTH', 'nominal', 1);
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, [area(3).both.pnom; area(1).both.pnom], 12, [t 'Pd']);
t_is(Qd, [area(3).both.qnom; area(1).both.qnom], 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone2, .type = ''BOTH'' .nominal = 0) : ';
opt = struct('type', 'BOTH', 'nominal', 0);
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, [area(3).both.p; area(1).both.p], 12, [t 'Pd']);
t_is(Qd, [area(3).both.q; area(1).both.q], 12, [t 'Qd']);

t = '      Pd = total_load(b, g, ld_zone2, .type = ''FIXED'') : ';
opt = struct('type', 'FIXED');
Pd = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, [area(3).fixed.p; area(1).fixed.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone2, .type = ''FIXED'') : ';
opt = struct('type', 'FIXED');
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, [area(3).fixed.p; area(1).fixed.p], 12, [t 'Pd']);
t_is(Qd, [area(3).fixed.q; area(1).fixed.q], 12, [t 'Qd']);

t = '      Pd = total_load(b, g, ld_zone2, .type = ''DISPATCHABLE'') : ';
opt = struct('type', 'DISPATCHABLE');
Pd = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, [area(3).disp.p; area(1).disp.p], 12, [t 'Pd']);

t = '      Pd = total_load(b, g, ld_zone2, .type = ''DISP'' .nominal = 1) : ';
opt = struct('type', 'DISP', 'nominal', 1);
Pd = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, [area(3).disp.pnom; area(1).disp.pnom], 12, [t 'Pd']);

t = '      Pd = total_load(b, g, ld_zone2, .type = ''DISP'' .nominal = 0) : ';
opt = struct('type', 'DISP', 'nominal', 0);
Pd = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, [area(3).disp.p; area(1).disp.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone2, .type = ''DISPATCHABLE'') : ';
opt = struct('type', 'DISPATCHABLE');
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, [area(3).disp.p; area(1).disp.p], 12, [t 'Pd']);
t_is(Qd, [area(3).disp.q; area(1).disp.q], 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone2, .type = ''DISP'' .nominal = 1) : ';
opt = struct('type', 'DISP', 'nominal', 1);
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, [area(3).disp.pnom; area(1).disp.pnom], 12, [t 'Pd']);
t_is(Qd, [area(3).disp.qnom; area(1).disp.qnom], 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone2, .type = ''DISP'' .nominal = 0) : ';
opt = struct('type', 'DISP', 'nominal', 0);
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, opt);
t_is(Pd, [area(3).disp.p; area(1).disp.p], 12, [t 'Pd']);
t_is(Qd, [area(3).disp.q; area(1).disp.q], 12, [t 'Qd']);

%%-----  old DEPRECATED string options, below  -----
%%-----  all load  -----
t = '      Pd = total_load(b, g, ''all'', ''BOTH'') : ';
Pd = total_load(mpc.bus, mpc.gen, 'all', 'BOTH');
t_is(Pd, total.both.pnom, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, ''all'', ''BOTH'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, 'all', 'BOTH');
t_is(Pd, total.both.pnom, 12, [t 'Pd']);
t_is(Qd, total.both.qnom, 12, [t 'Qd']);

t = '      Pd = total_load(b, g, ''all'', ''FIXED'') : ';
Pd = total_load(mpc.bus, mpc.gen, 'all', 'FIXED');
t_is(Pd, total.fixed.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, ''all'', ''FIXED'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, 'all', 'FIXED');
t_is(Pd, total.fixed.p, 12, [t 'Pd']);
t_is(Qd, total.fixed.q, 12, [t 'Qd']);

t = '      Pd = total_load(b, g, ''all'', ''DISPATCHABLE'') : ';
Pd = total_load(mpc.bus, mpc.gen, 'all', 'DISPATCHABLE');
t_is(Pd, total.disp.pnom, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, ''all'', ''DISPATCHABLE'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, 'all', 'DISPATCHABLE');
t_is(Pd, total.disp.pnom, 12, [t 'Pd']);
t_is(Qd, total.disp.qnom, 12, [t 'Qd']);

t = '      Pd = total_load(b, g, [], ''BOTH'') : ';
Pd = total_load(mpc.bus, mpc.gen, [], 'BOTH');
t_is(Pd, [area(1).both.pnom; area(2).both.pnom; area(3).both.pnom], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, [], ''BOTH'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, [], 'BOTH');
t_is(Pd, [area(1).both.pnom; area(2).both.pnom; area(3).both.pnom], 12, [t 'Pd']);
t_is(Qd, [area(1).both.qnom; area(2).both.qnom; area(3).both.qnom], 12, [t 'Qd']);

t = '      Pd = total_load(b, g, [], ''FIXED'') : ';
Pd = total_load(mpc.bus, mpc.gen, [], 'FIXED');
t_is(Pd, [area(1).fixed.p; area(2).fixed.p; area(3).fixed.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, [], ''FIXED'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, [], 'FIXED');
t_is(Pd, [area(1).fixed.p; area(2).fixed.p; area(3).fixed.p], 12, [t 'Pd']);
t_is(Qd, [area(1).fixed.q; area(2).fixed.q; area(3).fixed.q], 12, [t 'Qd']);

t = '      Pd = total_load(b, g, [], ''DISPATCHABLE'') : ';
Pd = total_load(mpc.bus, mpc.gen, [], 'DISPATCHABLE');
t_is(Pd, [area(1).disp.pnom; area(2).disp.pnom; area(3).disp.pnom], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, [], ''DISPATCHABLE'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, [], 'DISPATCHABLE');
t_is(Pd, [area(1).disp.pnom; area(2).disp.pnom; area(3).disp.pnom], 12, [t 'Pd']);
t_is(Qd, [area(1).disp.qnom; area(2).disp.qnom; area(3).disp.qnom], 12, [t 'Qd']);

%%-----  explicit single load zone  -----
nb = size(mpc.bus, 1);
load_zone = zeros(nb, 1);
k = find(mpc.bus(:, BUS_AREA) == 1);    %% area 1
load_zone(k) = 1;
t = '      Pd = total_load(b, g, ld_zone1, ''BOTH'') : ';
Pd = total_load(mpc.bus, mpc.gen, load_zone, 'BOTH');
t_is(Pd, area(1).both.pnom, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone1, ''BOTH'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, 'BOTH');
t_is(Pd, area(1).both.pnom, 12, [t 'Pd']);
t_is(Qd, area(1).both.qnom, 12, [t 'Qd']);

t = '      Pd = total_load(b, g, ld_zone1, ''FIXED'') : ';
Pd = total_load(mpc.bus, mpc.gen, load_zone, 'FIXED');
t_is(Pd, area(1).fixed.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone1, ''FIXED'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, 'FIXED');
t_is(Pd, area(1).fixed.p, 12, [t 'Pd']);
t_is(Qd, area(1).fixed.q, 12, [t 'Qd']);

t = '      Pd = total_load(b, g, ld_zone1, ''DISPATCHABLE'') : ';
Pd = total_load(mpc.bus, mpc.gen, load_zone, 'DISPATCHABLE');
t_is(Pd, area(1).disp.pnom, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone1, ''DISPATCHABLE'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, 'DISPATCHABLE');
t_is(Pd, area(1).disp.pnom, 12, [t 'Pd']);
t_is(Qd, area(1).disp.qnom, 12, [t 'Qd']);

%%-----  explicit multiple load zone  -----
load_zone = zeros(nb, 1);
k = find(mpc.bus(:, BUS_AREA) == 3);    %% area 3
load_zone(k) = 1;
k = find(mpc.bus(:, BUS_AREA) == 1);    %% area 1
load_zone(k) = 2;
t = '      Pd = total_load(b, g, ld_zone2, ''BOTH'') : ';
Pd = total_load(mpc.bus, mpc.gen, load_zone, 'BOTH');
t_is(Pd, [area(3).both.pnom; area(1).both.pnom], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone2, ''BOTH'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, 'BOTH');
t_is(Pd, [area(3).both.pnom; area(1).both.pnom], 12, [t 'Pd']);
t_is(Qd, [area(3).both.qnom; area(1).both.qnom], 12, [t 'Qd']);

t = '      Pd = total_load(b, g, ld_zone2, ''FIXED'') : ';
Pd = total_load(mpc.bus, mpc.gen, load_zone, 'FIXED');
t_is(Pd, [area(3).fixed.p; area(1).fixed.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone2, ''FIXED'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, 'FIXED');
t_is(Pd, [area(3).fixed.p; area(1).fixed.p], 12, [t 'Pd']);
t_is(Qd, [area(3).fixed.q; area(1).fixed.q], 12, [t 'Qd']);

t = '      Pd = total_load(b, g, ld_zone2, ''DISPATCHABLE'') : ';
Pd = total_load(mpc.bus, mpc.gen, load_zone, 'DISPATCHABLE');
t_is(Pd, [area(3).disp.pnom; area(1).disp.pnom], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(b, g, ld_zone2, ''DISPATCHABLE'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, 'DISPATCHABLE');
t_is(Pd, [area(3).disp.pnom; area(1).disp.pnom], 12, [t 'Pd']);
t_is(Qd, [area(3).disp.qnom; area(1).disp.qnom], 12, [t 'Qd']);

%%------------------------------------------------------------------
%% mostly same tests below, but with MPC input

 %%-----  all load  -----
t = '      Pd = total_load(mpc) : ';
Pd = total_load(mpc);
t_is(Pd, [area(1).both.p; area(2).both.p; area(3).both.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc) : ';
[Pd, Qd] = total_load(mpc);
t_is(Pd, [area(1).both.p; area(2).both.p; area(3).both.p], 12, [t 'Pd']);
t_is(Qd, [area(1).both.q; area(2).both.q; area(3).both.q], 12, [t 'Qd']);

t = '      Pd = total_load(mpc, ''all'') : ';
Pd = total_load(mpc, 'all');
t_is(Pd, total.both.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, ''all'') : ';
[Pd, Qd] = total_load(mpc, 'all');
t_is(Pd, total.both.p, 12, [t 'Pd']);
t_is(Qd, total.both.q, 12, [t 'Qd']);

t = '      Pd = total_load(mpc, ''all'', .type = ''BOTH'') : ';
opt = struct('type', 'BOTH');
Pd = total_load(mpc, 'all', opt);
t_is(Pd, total.both.p, 12, [t 'Pd']);

t = '      Pd = total_load(mpc, ''all'', .type = ''BOTH'' .nominal = 1) : ';
opt = struct('type', 'BOTH', 'nominal', 1);
Pd = total_load(mpc, 'all', opt);
t_is(Pd, total.both.pnom, 12, [t 'Pd']);

t = '      Pd = total_load(mpc, ''all'', .type = ''BOTH'' .nominal = 0) : ';
opt = struct('type', 'BOTH', 'nominal', 0);
Pd = total_load(mpc, 'all', opt);
t_is(Pd, total.both.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, ''all'', .type = ''BOTH'') : ';
opt = struct('type', 'BOTH');
[Pd, Qd] = total_load(mpc, 'all', opt);
t_is(Pd, total.both.p, 12, [t 'Pd']);
t_is(Qd, total.both.q, 12, [t 'Qd']);

t = '      Pd = total_load(mpc, ''all'', .type = ''FIXED'') : ';
opt = struct('type', 'FIXED');
Pd = total_load(mpc, 'all', opt);
t_is(Pd, total.fixed.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, ''all'', .type = ''FIXED'') : ';
opt = struct('type', 'FIXED');
[Pd, Qd] = total_load(mpc, 'all', opt);
t_is(Pd, total.fixed.p, 12, [t 'Pd']);
t_is(Qd, total.fixed.q, 12, [t 'Qd']);

t = '      Pd = total_load(mpc, ''all'', .type = ''DISPATCHABLE'') : ';
opt = struct('type', 'DISPATCHABLE');
Pd = total_load(mpc, 'all', opt);
t_is(Pd, total.disp.p, 12, [t 'Pd']);

t = '      Pd = total_load(mpc, ''all'', .type = ''DISP'' .nominal = 1) : ';
opt = struct('type', 'DISP', 'nominal', 1);
Pd = total_load(mpc, 'all', opt);
t_is(Pd, total.disp.pnom, 12, [t 'Pd']);

t = '      Pd = total_load(mpc, ''all'', .type = ''DISP'' .nominal = 0) : ';
opt = struct('type', 'DISP', 'nominal', 0);
Pd = total_load(mpc, 'all', opt);
t_is(Pd, total.disp.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, ''all'', .type = ''DISPATCHABLE'') : ';
opt = struct('type', 'DISPATCHABLE');
[Pd, Qd] = total_load(mpc, 'all', opt);
t_is(Pd, total.disp.p, 12, [t 'Pd']);
t_is(Qd, total.disp.q, 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(mpc, ''all'', .type = ''DISP'' .nominal = 1) : ';
opt = struct('type', 'DISP', 'nominal', 1);
[Pd, Qd] = total_load(mpc, 'all', opt);
t_is(Pd, total.disp.pnom, 12, [t 'Pd']);
t_is(Qd, total.disp.qnom, 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(mpc, ''all'', .type = ''DISP'' .nominal = 0) : ';
opt = struct('type', 'DISP', 'nominal', 0);
[Pd, Qd] = total_load(mpc, 'all', opt);
t_is(Pd, total.disp.p, 12, [t 'Pd']);
t_is(Qd, total.disp.q, 12, [t 'Qd']);

t = '      Pd = total_load(mpc, [], .type = ''BOTH'') : ';
opt = struct('type', 'BOTH');
Pd = total_load(mpc, [], opt);
t_is(Pd, [area(1).both.p; area(2).both.p; area(3).both.p], 12, [t 'Pd']);

t = '      Pd = total_load(mpc, [], .type = ''BOTH'' .nominal = 1) : ';
opt = struct('type', 'BOTH', 'nominal', 1);
Pd = total_load(mpc, [], opt);
t_is(Pd, [area(1).both.pnom; area(2).both.pnom; area(3).both.pnom], 12, [t 'Pd']);

t = '      Pd = total_load(mpc, [], .type = ''BOTH'' .nominal = 0) : ';
opt = struct('type', 'BOTH', 'nominal', 0);
Pd = total_load(mpc, [], opt);
t_is(Pd, [area(1).both.p; area(2).both.p; area(3).both.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, [], .type = ''BOTH'') : ';
opt = struct('type', 'BOTH');
[Pd, Qd] = total_load(mpc, [], opt);
t_is(Pd, [area(1).both.p; area(2).both.p; area(3).both.p], 12, [t 'Pd']);
t_is(Qd, [area(1).both.q; area(2).both.q; area(3).both.q], 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(mpc, [], .type = ''BOTH'' .nominal = 1) : ';
opt = struct('type', 'BOTH', 'nominal', 1);
[Pd, Qd] = total_load(mpc, [], opt);
t_is(Pd, [area(1).both.pnom; area(2).both.pnom; area(3).both.pnom], 12, [t 'Pd']);
t_is(Qd, [area(1).both.qnom; area(2).both.qnom; area(3).both.qnom], 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(mpc, [], .type = ''BOTH'' .nominal = 0) : ';
opt = struct('type', 'BOTH', 'nominal', 0);
[Pd, Qd] = total_load(mpc, [], opt);
t_is(Pd, [area(1).both.p; area(2).both.p; area(3).both.p], 12, [t 'Pd']);
t_is(Qd, [area(1).both.q; area(2).both.q; area(3).both.q], 12, [t 'Qd']);

t = '      Pd = total_load(mpc, [], .type = ''FIXED'') : ';
opt = struct('type', 'FIXED');
Pd = total_load(mpc, [], opt);
t_is(Pd, [area(1).fixed.p; area(2).fixed.p; area(3).fixed.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, [], .type = ''FIXED'') : ';
opt = struct('type', 'FIXED');
[Pd, Qd] = total_load(mpc, [], opt);
t_is(Pd, [area(1).fixed.p; area(2).fixed.p; area(3).fixed.p], 12, [t 'Pd']);
t_is(Qd, [area(1).fixed.q; area(2).fixed.q; area(3).fixed.q], 12, [t 'Qd']);

t = '      Pd = total_load(mpc, [], .type = ''DISPATCHABLE'') : ';
opt = struct('type', 'DISPATCHABLE');
Pd = total_load(mpc, [], opt);
t_is(Pd, [area(1).disp.p; area(2).disp.p; area(3).disp.p], 12, [t 'Pd']);

t = '      Pd = total_load(mpc, [], .type = ''DISP'' .nominal = 1) : ';
opt = struct('type', 'DISP', 'nominal', 1);
Pd = total_load(mpc, [], opt);
t_is(Pd, [area(1).disp.pnom; area(2).disp.pnom; area(3).disp.pnom], 12, [t 'Pd']);

t = '      Pd = total_load(mpc, [], .type = ''DISP'' .nominal = 0) : ';
opt = struct('type', 'DISP', 'nominal', 0);
Pd = total_load(mpc, [], opt);
t_is(Pd, [area(1).disp.p; area(2).disp.p; area(3).disp.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, [], .type = ''DISPATCHABLE'') : ';
opt = struct('type', 'DISPATCHABLE');
[Pd, Qd] = total_load(mpc, [], opt);
t_is(Pd, [area(1).disp.p; area(2).disp.p; area(3).disp.p], 12, [t 'Pd']);
t_is(Qd, [area(1).disp.q; area(2).disp.q; area(3).disp.q], 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(mpc, [], .type = ''DISP'' .nominal = 1) : ';
opt = struct('type', 'DISP', 'nominal', 1);
[Pd, Qd] = total_load(mpc, [], opt);
t_is(Pd, [area(1).disp.pnom; area(2).disp.pnom; area(3).disp.pnom], 12, [t 'Pd']);
t_is(Qd, [area(1).disp.qnom; area(2).disp.qnom; area(3).disp.qnom], 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(mpc, [], .type = ''DISP'' .nominal = 0) : ';
opt = struct('type', 'DISP', 'nominal', 0);
[Pd, Qd] = total_load(mpc, [], opt);
t_is(Pd, [area(1).disp.p; area(2).disp.p; area(3).disp.p], 12, [t 'Pd']);
t_is(Qd, [area(1).disp.q; area(2).disp.q; area(3).disp.q], 12, [t 'Qd']);

%%-----  explicit single load zone  -----
nb = size(mpc.bus, 1);
load_zone = zeros(nb, 1);
k = find(mpc.bus(:, BUS_AREA) == 1);    %% area 1
load_zone(k) = 1;
t = '      Pd = total_load(mpc, ld_zone1, .type = ''BOTH'') : ';
opt = struct('type', 'BOTH');
Pd = total_load(mpc, load_zone, opt);
t_is(Pd, area(1).both.p, 12, [t 'Pd']);

t = '      Pd = total_load(mpc, ld_zone1, .type = ''BOTH'' .nominal = 1) : ';
opt = struct('type', 'BOTH', 'nominal', 1);
Pd = total_load(mpc, load_zone, opt);
t_is(Pd, area(1).both.pnom, 12, [t 'Pd']);

t = '      Pd = total_load(mpc, ld_zone1, .type = ''BOTH'' .nominal = 0) : ';
opt = struct('type', 'BOTH', 'nominal', 0);
Pd = total_load(mpc, load_zone, opt);
t_is(Pd, area(1).both.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone1, .type = ''BOTH'') : ';
opt = struct('type', 'BOTH');
[Pd, Qd] = total_load(mpc, load_zone, opt);
t_is(Pd, area(1).both.p, 12, [t 'Pd']);
t_is(Qd, area(1).both.q, 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone1, .type = ''BOTH'' .nominal = 1) : ';
opt = struct('type', 'BOTH', 'nominal', 1);
[Pd, Qd] = total_load(mpc, load_zone, opt);
t_is(Pd, area(1).both.pnom, 12, [t 'Pd']);
t_is(Qd, area(1).both.qnom, 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone1, .type = ''BOTH'' .nominal = 0) : ';
opt = struct('type', 'BOTH', 'nominal', 0);
[Pd, Qd] = total_load(mpc, load_zone, opt);
t_is(Pd, area(1).both.p, 12, [t 'Pd']);
t_is(Qd, area(1).both.q, 12, [t 'Qd']);

t = '      Pd = total_load(mpc, ld_zone1, .type = ''FIXED'') : ';
opt = struct('type', 'FIXED');
Pd = total_load(mpc, load_zone, opt);
t_is(Pd, area(1).fixed.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone1, .type = ''FIXED'') : ';
opt = struct('type', 'FIXED');
[Pd, Qd] = total_load(mpc, load_zone, opt);
t_is(Pd, area(1).fixed.p, 12, [t 'Pd']);
t_is(Qd, area(1).fixed.q, 12, [t 'Qd']);

t = '      Pd = total_load(mpc, ld_zone1, .type = ''DISPATCHABLE'') : ';
opt = struct('type', 'DISPATCHABLE');
Pd = total_load(mpc, load_zone, opt);
t_is(Pd, area(1).disp.p, 12, [t 'Pd']);

t = '      Pd = total_load(mpc, ld_zone1, .type = ''DISP'' .nominal = 1) : ';
opt = struct('type', 'DISP', 'nominal', 1);
Pd = total_load(mpc, load_zone, opt);
t_is(Pd, area(1).disp.pnom, 12, [t 'Pd']);

t = '      Pd = total_load(mpc, ld_zone1, .type = ''DISP'' .nominal = 0) : ';
opt = struct('type', 'DISP', 'nominal', 0);
Pd = total_load(mpc, load_zone, opt);
t_is(Pd, area(1).disp.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone1, .type = ''DISPATCHABLE'') : ';
opt = struct('type', 'DISPATCHABLE');
[Pd, Qd] = total_load(mpc, load_zone, opt);
t_is(Pd, area(1).disp.p, 12, [t 'Pd']);
t_is(Qd, area(1).disp.q, 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone1, .type = ''DISP'' .nominal = 1) : ';
opt = struct('type', 'DISP', 'nominal', 1);
[Pd, Qd] = total_load(mpc, load_zone, opt);
t_is(Pd, area(1).disp.pnom, 12, [t 'Pd']);
t_is(Qd, area(1).disp.qnom, 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone1, .type = ''DISP'' .nominal = 0) : ';
opt = struct('type', 'DISP', 'nominal', 0);
[Pd, Qd] = total_load(mpc, load_zone, opt);
t_is(Pd, area(1).disp.p, 12, [t 'Pd']);
t_is(Qd, area(1).disp.q, 12, [t 'Qd']);

%%-----  explicit multiple load zone  -----
load_zone = zeros(nb, 1);
k = find(mpc.bus(:, BUS_AREA) == 3);    %% area 3
load_zone(k) = 1;
k = find(mpc.bus(:, BUS_AREA) == 1);    %% area 1
load_zone(k) = 2;
t = '      Pd = total_load(mpc, ld_zone2, .type = ''BOTH'') : ';
opt = struct('type', 'BOTH');
Pd = total_load(mpc, load_zone, opt);
t_is(Pd, [area(3).both.p; area(1).both.p], 12, [t 'Pd']);

t = '      Pd = total_load(mpc, ld_zone2, .type = ''BOTH'' .nominal = 1) : ';
opt = struct('type', 'BOTH', 'nominal', 1);
Pd = total_load(mpc, load_zone, opt);
t_is(Pd, [area(3).both.pnom; area(1).both.pnom], 12, [t 'Pd']);

t = '      Pd = total_load(mpc, ld_zone2, .type = ''BOTH'' .nominal = 0) : ';
opt = struct('type', 'BOTH', 'nominal', 0);
Pd = total_load(mpc, load_zone, opt);
t_is(Pd, [area(3).both.p; area(1).both.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone2, .type = ''BOTH'') : ';
opt = struct('type', 'BOTH');
[Pd, Qd] = total_load(mpc, load_zone, opt);
t_is(Pd, [area(3).both.p; area(1).both.p], 12, [t 'Pd']);
t_is(Qd, [area(3).both.q; area(1).both.q], 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone2, .type = ''BOTH'' .nominal = 1) : ';
opt = struct('type', 'BOTH', 'nominal', 1);
[Pd, Qd] = total_load(mpc, load_zone, opt);
t_is(Pd, [area(3).both.pnom; area(1).both.pnom], 12, [t 'Pd']);
t_is(Qd, [area(3).both.qnom; area(1).both.qnom], 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone2, .type = ''BOTH'' .nominal = 0) : ';
opt = struct('type', 'BOTH', 'nominal', 0);
[Pd, Qd] = total_load(mpc, load_zone, opt);
t_is(Pd, [area(3).both.p; area(1).both.p], 12, [t 'Pd']);
t_is(Qd, [area(3).both.q; area(1).both.q], 12, [t 'Qd']);

t = '      Pd = total_load(mpc, ld_zone2, .type = ''FIXED'') : ';
opt = struct('type', 'FIXED');
Pd = total_load(mpc, load_zone, opt);
t_is(Pd, [area(3).fixed.p; area(1).fixed.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone2, .type = ''FIXED'') : ';
opt = struct('type', 'FIXED');
[Pd, Qd] = total_load(mpc, load_zone, opt);
t_is(Pd, [area(3).fixed.p; area(1).fixed.p], 12, [t 'Pd']);
t_is(Qd, [area(3).fixed.q; area(1).fixed.q], 12, [t 'Qd']);

t = '      Pd = total_load(mpc, ld_zone2, .type = ''DISPATCHABLE'') : ';
opt = struct('type', 'DISPATCHABLE');
Pd = total_load(mpc, load_zone, opt);
t_is(Pd, [area(3).disp.p; area(1).disp.p], 12, [t 'Pd']);

t = '      Pd = total_load(mpc, ld_zone2, .type = ''DISP'' .nominal = 1) : ';
opt = struct('type', 'DISP', 'nominal', 1);
Pd = total_load(mpc, load_zone, opt);
t_is(Pd, [area(3).disp.pnom; area(1).disp.pnom], 12, [t 'Pd']);

t = '      Pd = total_load(mpc, ld_zone2, .type = ''DISP'' .nominal = 0) : ';
opt = struct('type', 'DISP', 'nominal', 0);
Pd = total_load(mpc, load_zone, opt);
t_is(Pd, [area(3).disp.p; area(1).disp.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone2, .type = ''DISPATCHABLE'') : ';
opt = struct('type', 'DISPATCHABLE');
[Pd, Qd] = total_load(mpc, load_zone, opt);
t_is(Pd, [area(3).disp.p; area(1).disp.p], 12, [t 'Pd']);
t_is(Qd, [area(3).disp.q; area(1).disp.q], 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone2, .type = ''DISP'' .nominal = 1) : ';
opt = struct('type', 'DISP', 'nominal', 1);
[Pd, Qd] = total_load(mpc, load_zone, opt);
t_is(Pd, [area(3).disp.pnom; area(1).disp.pnom], 12, [t 'Pd']);
t_is(Qd, [area(3).disp.qnom; area(1).disp.qnom], 12, [t 'Qd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone2, .type = ''DISP'' .nominal = 0) : ';
opt = struct('type', 'DISP', 'nominal', 0);
[Pd, Qd] = total_load(mpc, load_zone, opt);
t_is(Pd, [area(3).disp.p; area(1).disp.p], 12, [t 'Pd']);
t_is(Qd, [area(3).disp.q; area(1).disp.q], 12, [t 'Qd']);

%%-----  old DEPRECATED string options, below  -----
%%-----  all load  -----
t = '      Pd = total_load(mpc, ''all'', ''BOTH'') : ';
Pd = total_load(mpc, 'all', 'BOTH');
t_is(Pd, total.both.pnom, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, ''all'', ''BOTH'') : ';
[Pd, Qd] = total_load(mpc, 'all', 'BOTH');
t_is(Pd, total.both.pnom, 12, [t 'Pd']);
t_is(Qd, total.both.qnom, 12, [t 'Qd']);

t = '      Pd = total_load(mpc, ''all'', ''FIXED'') : ';
Pd = total_load(mpc, 'all', 'FIXED');
t_is(Pd, total.fixed.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, ''all'', ''FIXED'') : ';
[Pd, Qd] = total_load(mpc, 'all', 'FIXED');
t_is(Pd, total.fixed.p, 12, [t 'Pd']);
t_is(Qd, total.fixed.q, 12, [t 'Qd']);

t = '      Pd = total_load(mpc, ''all'', ''DISPATCHABLE'') : ';
Pd = total_load(mpc, 'all', 'DISPATCHABLE');
t_is(Pd, total.disp.pnom, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, ''all'', ''DISPATCHABLE'') : ';
[Pd, Qd] = total_load(mpc, 'all', 'DISPATCHABLE');
t_is(Pd, total.disp.pnom, 12, [t 'Pd']);
t_is(Qd, total.disp.qnom, 12, [t 'Qd']);

t = '      Pd = total_load(mpc, [], ''BOTH'') : ';
Pd = total_load(mpc, [], 'BOTH');
t_is(Pd, [area(1).both.pnom; area(2).both.pnom; area(3).both.pnom], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, [], ''BOTH'') : ';
[Pd, Qd] = total_load(mpc, [], 'BOTH');
t_is(Pd, [area(1).both.pnom; area(2).both.pnom; area(3).both.pnom], 12, [t 'Pd']);
t_is(Qd, [area(1).both.qnom; area(2).both.qnom; area(3).both.qnom], 12, [t 'Qd']);

t = '      Pd = total_load(mpc, [], ''FIXED'') : ';
Pd = total_load(mpc, [], 'FIXED');
t_is(Pd, [area(1).fixed.p; area(2).fixed.p; area(3).fixed.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, [], ''FIXED'') : ';
[Pd, Qd] = total_load(mpc, [], 'FIXED');
t_is(Pd, [area(1).fixed.p; area(2).fixed.p; area(3).fixed.p], 12, [t 'Pd']);
t_is(Qd, [area(1).fixed.q; area(2).fixed.q; area(3).fixed.q], 12, [t 'Qd']);

t = '      Pd = total_load(mpc, [], ''DISPATCHABLE'') : ';
Pd = total_load(mpc, [], 'DISPATCHABLE');
t_is(Pd, [area(1).disp.pnom; area(2).disp.pnom; area(3).disp.pnom], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, [], ''DISPATCHABLE'') : ';
[Pd, Qd] = total_load(mpc, [], 'DISPATCHABLE');
t_is(Pd, [area(1).disp.pnom; area(2).disp.pnom; area(3).disp.pnom], 12, [t 'Pd']);
t_is(Qd, [area(1).disp.qnom; area(2).disp.qnom; area(3).disp.qnom], 12, [t 'Qd']);

%%-----  explicit single load zone  -----
nb = size(mpc.bus, 1);
load_zone = zeros(nb, 1);
k = find(mpc.bus(:, BUS_AREA) == 1);    %% area 1
load_zone(k) = 1;
t = '      Pd = total_load(mpc, ld_zone1, ''BOTH'') : ';
Pd = total_load(mpc, load_zone, 'BOTH');
t_is(Pd, area(1).both.pnom, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone1, ''BOTH'') : ';
[Pd, Qd] = total_load(mpc, load_zone, 'BOTH');
t_is(Pd, area(1).both.pnom, 12, [t 'Pd']);
t_is(Qd, area(1).both.qnom, 12, [t 'Qd']);

t = '      Pd = total_load(mpc, ld_zone1, ''FIXED'') : ';
Pd = total_load(mpc, load_zone, 'FIXED');
t_is(Pd, area(1).fixed.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone1, ''FIXED'') : ';
[Pd, Qd] = total_load(mpc, load_zone, 'FIXED');
t_is(Pd, area(1).fixed.p, 12, [t 'Pd']);
t_is(Qd, area(1).fixed.q, 12, [t 'Qd']);

t = '      Pd = total_load(mpc, ld_zone1, ''DISPATCHABLE'') : ';
Pd = total_load(mpc, load_zone, 'DISPATCHABLE');
t_is(Pd, area(1).disp.pnom, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone1, ''DISPATCHABLE'') : ';
[Pd, Qd] = total_load(mpc, load_zone, 'DISPATCHABLE');
t_is(Pd, area(1).disp.pnom, 12, [t 'Pd']);
t_is(Qd, area(1).disp.qnom, 12, [t 'Qd']);

%%-----  explicit multiple load zone  -----
load_zone = zeros(nb, 1);
k = find(mpc.bus(:, BUS_AREA) == 3);    %% area 3
load_zone(k) = 1;
k = find(mpc.bus(:, BUS_AREA) == 1);    %% area 1
load_zone(k) = 2;
t = '      Pd = total_load(mpc, ld_zone2, ''BOTH'') : ';
Pd = total_load(mpc, load_zone, 'BOTH');
t_is(Pd, [area(3).both.pnom; area(1).both.pnom], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone2, ''BOTH'') : ';
[Pd, Qd] = total_load(mpc, load_zone, 'BOTH');
t_is(Pd, [area(3).both.pnom; area(1).both.pnom], 12, [t 'Pd']);
t_is(Qd, [area(3).both.qnom; area(1).both.qnom], 12, [t 'Qd']);

t = '      Pd = total_load(mpc, ld_zone2, ''FIXED'') : ';
Pd = total_load(mpc, load_zone, 'FIXED');
t_is(Pd, [area(3).fixed.p; area(1).fixed.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone2, ''FIXED'') : ';
[Pd, Qd] = total_load(mpc, load_zone, 'FIXED');
t_is(Pd, [area(3).fixed.p; area(1).fixed.p], 12, [t 'Pd']);
t_is(Qd, [area(3).fixed.q; area(1).fixed.q], 12, [t 'Qd']);

t = '      Pd = total_load(mpc, ld_zone2, ''DISPATCHABLE'') : ';
Pd = total_load(mpc, load_zone, 'DISPATCHABLE');
t_is(Pd, [area(3).disp.pnom; area(1).disp.pnom], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(mpc, ld_zone2, ''DISPATCHABLE'') : ';
[Pd, Qd] = total_load(mpc, load_zone, 'DISPATCHABLE');
t_is(Pd, [area(3).disp.pnom; area(1).disp.pnom], 12, [t 'Pd']);
t_is(Qd, [area(3).disp.qnom; area(1).disp.qnom], 12, [t 'Qd']);

t_end;
