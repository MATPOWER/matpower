function t_total_load(quiet)
%T_TOTAL_LOAD  Tests for code in TOTAL_LOAD.

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
%   other modules (such as M-files and MEX-files) available in a
%   Matlab (or compatible) environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

if nargin < 1
    quiet = 0;
end

n_tests = 48;

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

t = '      Pd = total_load(bus, gen, ''all'') : ';
Pd = total_load(mpc.bus, mpc.gen, 'all');
t_is(Pd, total.both.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(bus, gen, ''all'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, 'all');
t_is(Pd, total.both.p, 12, [t 'Pd']);
t_is(Qd, total.both.q, 12, [t 'Qd']);

t = '      Pd = total_load(bus, gen, ''all'', ''BOTH'') : ';
Pd = total_load(mpc.bus, mpc.gen, 'all', 'BOTH');
t_is(Pd, total.both.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(bus, gen, ''all'', ''BOTH'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, 'all', 'BOTH');
t_is(Pd, total.both.p, 12, [t 'Pd']);
t_is(Qd, total.both.q, 12, [t 'Qd']);

t = '      Pd = total_load(bus, gen, ''all'', ''FIXED'') : ';
Pd = total_load(mpc.bus, mpc.gen, 'all', 'FIXED');
t_is(Pd, total.fixed.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(bus, gen, ''all'', ''FIXED'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, 'all', 'FIXED');
t_is(Pd, total.fixed.p, 12, [t 'Pd']);
t_is(Qd, total.fixed.q, 12, [t 'Qd']);

t = '      Pd = total_load(bus, gen, ''all'', ''DISPATCHABLE'') : ';
Pd = total_load(mpc.bus, mpc.gen, 'all', 'DISPATCHABLE');
t_is(Pd, total.disp.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(bus, gen, ''all'', ''DISPATCHABLE'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, 'all', 'DISPATCHABLE');
t_is(Pd, total.disp.p, 12, [t 'Pd']);
t_is(Qd, total.disp.q, 12, [t 'Qd']);

t = '      Pd = total_load(bus, gen, [], ''BOTH'') : ';
Pd = total_load(mpc.bus, mpc.gen, [], 'BOTH');
t_is(Pd, [area(1).both.p; area(2).both.p; area(3).both.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(bus, gen, [], ''BOTH'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, [], 'BOTH');
t_is(Pd, [area(1).both.p; area(2).both.p; area(3).both.p], 12, [t 'Pd']);
t_is(Qd, [area(1).both.q; area(2).both.q; area(3).both.q], 12, [t 'Qd']);

t = '      Pd = total_load(bus, gen, [], ''FIXED'') : ';
Pd = total_load(mpc.bus, mpc.gen, [], 'FIXED');
t_is(Pd, [area(1).fixed.p; area(2).fixed.p; area(3).fixed.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(bus, gen, [], ''FIXED'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, [], 'FIXED');
t_is(Pd, [area(1).fixed.p; area(2).fixed.p; area(3).fixed.p], 12, [t 'Pd']);
t_is(Qd, [area(1).fixed.q; area(2).fixed.q; area(3).fixed.q], 12, [t 'Qd']);

t = '      Pd = total_load(bus, gen, [], ''DISPATCHABLE'') : ';
Pd = total_load(mpc.bus, mpc.gen, [], 'DISPATCHABLE');
t_is(Pd, [area(1).disp.p; area(2).disp.p; area(3).disp.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(bus, gen, [], ''DISPATCHABLE'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, [], 'DISPATCHABLE');
t_is(Pd, [area(1).disp.p; area(2).disp.p; area(3).disp.p], 12, [t 'Pd']);
t_is(Qd, [area(1).disp.q; area(2).disp.q; area(3).disp.q], 12, [t 'Qd']);

%%-----  explicit single load zone  -----
nb = size(mpc.bus, 1);
load_zone = zeros(nb, 1);
k = find(mpc.bus(:, BUS_AREA) == 2);    %% area 2
load_zone(k) = 1;
t = '      Pd = total_load(bus, gen, load_zone1, ''BOTH'') : ';
Pd = total_load(mpc.bus, mpc.gen, load_zone, 'BOTH');
t_is(Pd, area(2).both.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(bus, gen, load_zone1, ''BOTH'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, 'BOTH');
t_is(Pd, area(2).both.p, 12, [t 'Pd']);
t_is(Qd, area(2).both.q, 12, [t 'Qd']);

t = '      Pd = total_load(bus, gen, load_zone1, ''FIXED'') : ';
Pd = total_load(mpc.bus, mpc.gen, load_zone, 'FIXED');
t_is(Pd, area(2).fixed.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(bus, gen, load_zone1, ''FIXED'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, 'FIXED');
t_is(Pd, area(2).fixed.p, 12, [t 'Pd']);
t_is(Qd, area(2).fixed.q, 12, [t 'Qd']);

t = '      Pd = total_load(bus, gen, load_zone1, ''DISPATCHABLE'') : ';
Pd = total_load(mpc.bus, mpc.gen, load_zone, 'DISPATCHABLE');
t_is(Pd, area(2).disp.p, 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(bus, gen, load_zone1, ''DISPATCHABLE'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, 'DISPATCHABLE');
t_is(Pd, area(2).disp.p, 12, [t 'Pd']);
t_is(Qd, area(2).disp.q, 12, [t 'Qd']);

%%-----  explicit multiple load zone  -----
load_zone = zeros(nb, 1);
k = find(mpc.bus(:, BUS_AREA) == 3);    %% area 3
load_zone(k) = 1;
k = find(mpc.bus(:, BUS_AREA) == 1);    %% area 1
load_zone(k) = 2;
t = '      Pd = total_load(bus, gen, load_zone2, ''BOTH'') : ';
Pd = total_load(mpc.bus, mpc.gen, load_zone, 'BOTH');
t_is(Pd, [area(3).both.p; area(1).both.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(bus, gen, load_zone2, ''BOTH'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, 'BOTH');
t_is(Pd, [area(3).both.p; area(1).both.p], 12, [t 'Pd']);
t_is(Qd, [area(3).both.q; area(1).both.q], 12, [t 'Qd']);

t = '      Pd = total_load(bus, gen, load_zone2, ''FIXED'') : ';
Pd = total_load(mpc.bus, mpc.gen, load_zone, 'FIXED');
t_is(Pd, [area(3).fixed.p; area(1).fixed.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(bus, gen, load_zone2, ''FIXED'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, 'FIXED');
t_is(Pd, [area(3).fixed.p; area(1).fixed.p], 12, [t 'Pd']);
t_is(Qd, [area(3).fixed.q; area(1).fixed.q], 12, [t 'Qd']);

t = '      Pd = total_load(bus, gen, load_zone2, ''DISPATCHABLE'') : ';
Pd = total_load(mpc.bus, mpc.gen, load_zone, 'DISPATCHABLE');
t_is(Pd, [area(3).disp.p; area(1).disp.p], 12, [t 'Pd']);

t = '[Pd, Qd] = total_load(bus, gen, load_zone2, ''DISPATCHABLE'') : ';
[Pd, Qd] = total_load(mpc.bus, mpc.gen, load_zone, 'DISPATCHABLE');
t_is(Pd, [area(3).disp.p; area(1).disp.p], 12, [t 'Pd']);
t_is(Qd, [area(3).disp.q; area(1).disp.q], 12, [t 'Qd']);

t_end;
