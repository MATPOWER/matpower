function t_apply_changes(quiet)
%T_APPLY_CHANGES  Tests for code in apply_changes.m.

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

n_tests = 570;

t_begin(n_tests, quiet);

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
[CT_LABEL, CT_PROB, CT_TABLE, CT_TBUS, CT_TGEN, CT_TBRCH, CT_TAREABUS, ...
    CT_TAREAGEN, CT_TAREABRCH, CT_ROW, CT_COL, CT_CHGTYPE, CT_REP, ...
    CT_REL, CT_ADD, CT_NEWVAL, CT_TLOAD, CT_TAREALOAD, CT_LOAD_ALL_PQ, ...
    CT_LOAD_FIX_PQ, CT_LOAD_DIS_PQ, CT_LOAD_ALL_P, CT_LOAD_FIX_P, ...
    CT_LOAD_DIS_P, CT_TGENCOST, CT_TAREAGENCOST, CT_MODCOST_F, ...
    CT_MODCOST_X] = idx_ct;

mpc = loadcase('t_auction_case');
mpc.gen(8, GEN_BUS) = 2;    %% multiple d. loads per area, same bus as gen
mpc.gen(8, [QG QMIN QMAX]) = [ 3 0 3 ];
%% put it load before gen in matrix
mpc.gen = [mpc.gen(8, :); mpc.gen(1:7, :); mpc.gen(9, :)];
ld = find(isload(mpc.gen));

%% create map of external bus numbers to bus indices
i2e = mpc.bus(:, BUS_I);
e2i = sparse(max(i2e), 1);
e2i(i2e) = (1:size(mpc.bus, 1))';

for k = 1:3
    a{k} = find(mpc.bus(:, BUS_AREA) == k); %% buses in area k
    [junk, tmp, junk2] = intersect(e2i(mpc.gen(ld, GEN_BUS)), a{k});
    lda{k} = ld(tmp);                       %% disp loads in area k
    ga{k} = find(ismember(e2i(mpc.gen(:, GEN_BUS)), a{k}));
    tmp1 = find(ismember(e2i(mpc.branch(:, F_BUS)), a{k}));
    tmp2 = find(ismember(e2i(mpc.branch(:, T_BUS)), a{k}));
    bra{k} = unique([tmp1; tmp2]);
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

chgtab = [
    1   0.002   CT_TBUS         1   GS              CT_REP  5;
    2   0.002   CT_TBUS         2   VMIN            CT_REL  1/0.95; %% 1
    3   0.002   CT_TBUS         0   VMAX            CT_ADD  -0.01;  %% VMAX - 0.01
    4   0.002   CT_TBRCH        1   BR_STATUS       CT_REP  0;
    5   0.002   CT_TBRCH        2   RATE_A          CT_REL  0.1;    %% 1000
    6   0.002   CT_TBRCH        0   RATE_B          CT_ADD  0.1;    %% RATE_B + 0.1
    7   0.002   CT_TGEN         1   GEN_STATUS      CT_REP  0;
    8   0.002   CT_TGEN         2   QMAX            CT_REL  1.1;    %% 66
    9   0.002   CT_TGEN         0   PMIN            CT_ADD  0.5;    %% PMIN + 0.5
    10  0.002   CT_TAREABUS     1   VMAX            CT_REP  1.1;
    11  0.002   CT_TAREABUS     2   VMAX            CT_REL  1.01;
    12  0.002   CT_TAREABUS     3   VMAX            CT_ADD  0.1;
    13  0.002   CT_TAREABRCH    1   RATE_B          CT_REP  100;
    14  0.002   CT_TAREABRCH    2   RATE_B          CT_REL  1.1;
    15  0.002   CT_TAREABRCH    3   RATE_B          CT_ADD  0.5;
    16  0.002   CT_TAREAGEN     1   PMIN            CT_REP  0;
    17  0.002   CT_TAREAGEN     2   PMIN            CT_REL  1.1;
    18  0.002   CT_TAREAGEN     3   PMIN            CT_ADD  0.1;
    20  0.002   CT_TLOAD        0   CT_LOAD_FIX_PQ  CT_REL  2;      %% fixed PQ * 2
    21  0.002   CT_TLOAD        0   CT_LOAD_FIX_P   CT_REL  2;      %% fixed P * 2
    22  0.002   CT_TLOAD        0   CT_LOAD_ALL_PQ  CT_REL  2;      %% all PQ * 2
    23  0.002   CT_TLOAD        0   CT_LOAD_ALL_P   CT_REL  2;      %% all P * 2
    24  0.002   CT_TLOAD        0   CT_LOAD_DIS_PQ  CT_REL  2;      %% disp PQ * 2
    25  0.002   CT_TLOAD        0   CT_LOAD_DIS_P   CT_REL  2;      %% disp P * 2
    26  0.002   CT_TLOAD        0   CT_LOAD_FIX_PQ  CT_REP  200;    %% fixed PQ => total = 200
    27  0.002   CT_TLOAD        0   CT_LOAD_FIX_P   CT_REP  200;    %% fixed P => total = 200
    28  0.002   CT_TLOAD        0   CT_LOAD_ALL_PQ  CT_REP  200;    %% all PQ => total = 200
    29  0.002   CT_TLOAD        0   CT_LOAD_ALL_P   CT_REP  200;    %% all P => total = 200
    30  0.002   CT_TLOAD        0   CT_LOAD_DIS_PQ  CT_REP  200;    %% disp PQ => total = 200
    31  0.002   CT_TLOAD        0   CT_LOAD_DIS_P   CT_REP  200;    %% disp P => total = 200
    32  0.002   CT_TLOAD        0   CT_LOAD_FIX_PQ  CT_ADD  25;     %% fixed PQ + 25
    33  0.002   CT_TLOAD        0   CT_LOAD_FIX_P   CT_ADD  25;     %% fixed P + 25
    34  0.002   CT_TLOAD        0   CT_LOAD_ALL_PQ  CT_ADD  25;     %% all PQ + 25
    35  0.002   CT_TLOAD        0   CT_LOAD_ALL_P   CT_ADD  25;     %% all P + 25
    36  0.002   CT_TLOAD        0   CT_LOAD_DIS_PQ  CT_ADD  25;     %% disp PQ + 25
    37  0.002   CT_TLOAD        0   CT_LOAD_DIS_P   CT_ADD  25;     %% disp P + 25
    40  0.002   CT_TLOAD        2   CT_LOAD_FIX_PQ  CT_REL  2;      %% fixed PQ * 2
    41  0.002   CT_TLOAD        2   CT_LOAD_FIX_P   CT_REL  2;      %% fixed P * 2
    42  0.002   CT_TLOAD        2   CT_LOAD_ALL_PQ  CT_REL  2;      %% all PQ * 2
    43  0.002   CT_TLOAD        2   CT_LOAD_ALL_P   CT_REL  2;      %% all P * 2
    44  0.002   CT_TLOAD        2   CT_LOAD_DIS_PQ  CT_REL  2;      %% disp PQ * 2
    45  0.002   CT_TLOAD        2   CT_LOAD_DIS_P   CT_REL  2;      %% disp P * 2
    46  0.002   CT_TLOAD        2   CT_LOAD_FIX_PQ  CT_REP  50;     %% fixed PQ => total = 50
    47  0.002   CT_TLOAD        2   CT_LOAD_FIX_P   CT_REP  50;     %% fixed P => total = 50
    48  0.002   CT_TLOAD        2   CT_LOAD_ALL_PQ  CT_REP  50;     %% all PQ => total = 50
    49  0.002   CT_TLOAD        2   CT_LOAD_ALL_P   CT_REP  50;     %% all P => total = 50
    50  0.002   CT_TLOAD        2   CT_LOAD_DIS_PQ  CT_REP  50;     %% disp PQ => total = 50
    51  0.002   CT_TLOAD        2   CT_LOAD_DIS_P   CT_REP  50;     %% disp P => total = 50
    52  0.002   CT_TLOAD        2   CT_LOAD_FIX_PQ  CT_ADD  10;     %% fixed PQ + 10
    53  0.002   CT_TLOAD        2   CT_LOAD_FIX_P   CT_ADD  10;     %% fixed P + 10
    54  0.002   CT_TLOAD        2   CT_LOAD_ALL_PQ  CT_ADD  10;     %% all PQ + 10
    55  0.002   CT_TLOAD        2   CT_LOAD_ALL_P   CT_ADD  10;     %% all P + 10
    56  0.002   CT_TLOAD        2   CT_LOAD_DIS_PQ  CT_ADD  10;     %% disp PQ + 10
    57  0.002   CT_TLOAD        2   CT_LOAD_DIS_P   CT_ADD  10;     %% disp P + 10
    60  0.002   CT_TAREALOAD    1   CT_LOAD_FIX_PQ  CT_REL  3;      %% fixed PQ * 2
    61  0.002   CT_TAREALOAD    1   CT_LOAD_FIX_P   CT_REL  3;      %% fixed P * 2
    62  0.002   CT_TAREALOAD    1   CT_LOAD_ALL_PQ  CT_REL  3;      %% all PQ * 2
    63  0.002   CT_TAREALOAD    1   CT_LOAD_ALL_P   CT_REL  3;      %% all P * 2
    64  0.002   CT_TAREALOAD    1   CT_LOAD_DIS_PQ  CT_REL  3;      %% disp PQ * 2
    65  0.002   CT_TAREALOAD    1   CT_LOAD_DIS_P   CT_REL  3;      %% disp P * 2
    66  0.002   CT_TAREALOAD    1   CT_LOAD_FIX_PQ  CT_REP  100;    %% fixed PQ => total = 50
    67  0.002   CT_TAREALOAD    1   CT_LOAD_FIX_P   CT_REP  100;    %% fixed P => total = 50
    68  0.002   CT_TAREALOAD    1   CT_LOAD_ALL_PQ  CT_REP  100;    %% all PQ => total = 50
    69  0.002   CT_TAREALOAD    1   CT_LOAD_ALL_P   CT_REP  100;    %% all P => total = 50
    70  0.002   CT_TAREALOAD    1   CT_LOAD_DIS_PQ  CT_REP  100;    %% disp PQ => total = 50
    71  0.002   CT_TAREALOAD    1   CT_LOAD_DIS_P   CT_REP  100;    %% disp P => total = 50
    72  0.002   CT_TAREALOAD    1   CT_LOAD_FIX_PQ  CT_ADD  20;     %% fixed PQ + 10
    73  0.002   CT_TAREALOAD    1   CT_LOAD_FIX_P   CT_ADD  20;     %% fixed P + 10
    74  0.002   CT_TAREALOAD    1   CT_LOAD_ALL_PQ  CT_ADD  20;     %% all PQ + 10
    75  0.002   CT_TAREALOAD    1   CT_LOAD_ALL_P   CT_ADD  20;     %% all P + 10
    76  0.002   CT_TAREALOAD    1   CT_LOAD_DIS_PQ  CT_ADD  20;     %% disp PQ + 10
    77  0.002   CT_TAREALOAD    1   CT_LOAD_DIS_P   CT_ADD  20;     %% disp P + 10
    78  0.002   CT_TGENCOST     1   STARTUP         CT_REP  1000;
    79  0.002   CT_TGENCOST     2   COST+3          CT_REL  1.1;    %% 240 * 1.1
    80  0.002   CT_TGENCOST     0   NCOST           CT_ADD  -1;     %% NCOST - 1
    81  0.002   CT_TGENCOST     1   CT_MODCOST_X    CT_REL  1.1;    %% SCALE COST HORIZ
    82  0.002   CT_TGENCOST     0   CT_MODCOST_X    CT_REL  0.9;    %% SCALE COST HORIZ
    83  0.002   CT_TGENCOST     2   CT_MODCOST_X    CT_ADD  -10;    %% SHIFT COST HORIZ
    84  0.002   CT_TGENCOST     0   CT_MODCOST_X    CT_ADD  20;     %% SHIFT COST HORIZ
    85  0.002   CT_TGENCOST     1   CT_MODCOST_F    CT_REL  1.1;    %% SCALE COST VERTICAL
    86  0.002   CT_TGENCOST     0   CT_MODCOST_F    CT_REL  0.9;    %% SCALE COST VERTICAL
    87  0.002   CT_TGENCOST     2   CT_MODCOST_F    CT_ADD  -10;    %% SHIFT COST VERTICAL
    88  0.002   CT_TGENCOST     0   CT_MODCOST_F    CT_ADD  20;     %% SHIFT COST VERTICAL
    89  0.002   CT_TAREAGENCOST 1   SHUTDOWN        CT_REP  500;
    90  0.002   CT_TAREAGENCOST 2   COST+5          CT_REL  1.1;
    91  0.002   CT_TAREAGENCOST 3   NCOST           CT_ADD  -2;
    92  0.002   CT_TAREAGENCOST 1   CT_MODCOST_X    CT_REL  1.1;    %% SCALE COST HORIZ
    93  0.002   CT_TAREAGENCOST 2   CT_MODCOST_X    CT_ADD  -10;    %% SHIFT COST HORIZ
    94  0.002   CT_TAREAGENCOST 1   CT_MODCOST_F    CT_REL  1.1;    %% SCALE COST VERTICAL
    95  0.002   CT_TAREAGENCOST 2   CT_MODCOST_F    CT_ADD  -10;    %% SHIFT COST VERTICAL
    122 0.002   CT_TLOAD        0   -CT_LOAD_ALL_PQ CT_REL  2;      %% all PQ * 2, incl cost
    123 0.002   CT_TLOAD        0   -CT_LOAD_ALL_P  CT_REL  2;      %% all P * 2, incl cost
    124 0.002   CT_TLOAD        0   -CT_LOAD_DIS_PQ CT_REL  2;      %% disp PQ * 2, incl cost
    125 0.002   CT_TLOAD        0   -CT_LOAD_DIS_P  CT_REL  2;      %% disp P * 2, incl cost
    128 0.002   CT_TLOAD        0   -CT_LOAD_ALL_PQ CT_REP  200;    %% all PQ => total = 200, incl cost
    129 0.002   CT_TLOAD        0   -CT_LOAD_ALL_P  CT_REP  200;    %% all P => total = 200, incl cost
    130 0.002   CT_TLOAD        0   -CT_LOAD_DIS_PQ CT_REP  200;    %% disp PQ => total = 200, incl cost
    131 0.002   CT_TLOAD        0   -CT_LOAD_DIS_P  CT_REP  200;    %% disp P => total = 200, incl cost
    134 0.002   CT_TLOAD        0   -CT_LOAD_ALL_PQ CT_ADD  25;     %% all PQ + 25, incl cost
    135 0.002   CT_TLOAD        0   -CT_LOAD_ALL_P  CT_ADD  25;     %% all P + 25, incl cost
    136 0.002   CT_TLOAD        0   -CT_LOAD_DIS_PQ CT_ADD  25;     %% disp PQ + 25, incl cost
    137 0.002   CT_TLOAD        0   -CT_LOAD_DIS_P  CT_ADD  25;     %% disp P + 25, incl cost
    142 0.002   CT_TLOAD        2   -CT_LOAD_ALL_PQ CT_REL  2;      %% all PQ * 2, incl cost
    143 0.002   CT_TLOAD        2   -CT_LOAD_ALL_P  CT_REL  2;      %% all P * 2, incl cost
    144 0.002   CT_TLOAD        2   -CT_LOAD_DIS_PQ CT_REL  2;      %% disp PQ * 2, incl cost
    145 0.002   CT_TLOAD        2   -CT_LOAD_DIS_P  CT_REL  2;      %% disp P * 2, incl cost
    148 0.002   CT_TLOAD        2   -CT_LOAD_ALL_PQ CT_REP  50;     %% all PQ => total = 50, incl cost
    149 0.002   CT_TLOAD        2   -CT_LOAD_ALL_P  CT_REP  50;     %% all P => total = 50, incl cost
    150 0.002   CT_TLOAD        2   -CT_LOAD_DIS_PQ CT_REP  50;     %% disp PQ => total = 50, incl cost
    151 0.002   CT_TLOAD        2   -CT_LOAD_DIS_P  CT_REP  50;     %% disp P => total = 50, incl cost
    154 0.002   CT_TLOAD        2   -CT_LOAD_ALL_PQ CT_ADD  10;     %% all PQ + 10, incl cost
    155 0.002   CT_TLOAD        2   -CT_LOAD_ALL_P  CT_ADD  10;     %% all P + 10, incl cost
    156 0.002   CT_TLOAD        2   -CT_LOAD_DIS_PQ CT_ADD  10;     %% disp PQ + 10, incl cost
    157 0.002   CT_TLOAD        2   -CT_LOAD_DIS_P  CT_ADD  10;     %% disp P + 10, incl cost
    162 0.002   CT_TAREALOAD    1   -CT_LOAD_ALL_PQ CT_REL  3;      %% all PQ * 2, incl cost
    163 0.002   CT_TAREALOAD    1   -CT_LOAD_ALL_P  CT_REL  3;      %% all P * 2, incl cost
    164 0.002   CT_TAREALOAD    1   -CT_LOAD_DIS_PQ CT_REL  3;      %% disp PQ * 2, incl cost
    165 0.002   CT_TAREALOAD    1   -CT_LOAD_DIS_P  CT_REL  3;      %% disp P * 2, incl cost
    168 0.002   CT_TAREALOAD    1   -CT_LOAD_ALL_PQ CT_REP  100;    %% all PQ => total = 50, incl cost
    169 0.002   CT_TAREALOAD    1   -CT_LOAD_ALL_P  CT_REP  100;    %% all P => total = 50, incl cost
    170 0.002   CT_TAREALOAD    1   -CT_LOAD_DIS_PQ CT_REP  100;    %% disp PQ => total = 50, incl cost
    171 0.002   CT_TAREALOAD    1   -CT_LOAD_DIS_P  CT_REP  100;    %% disp P => total = 50, incl cost
    174 0.002   CT_TAREALOAD    1   -CT_LOAD_ALL_PQ CT_ADD  20;     %% all PQ + 10, incl cost
    175 0.002   CT_TAREALOAD    1   -CT_LOAD_ALL_P  CT_ADD  20;     %% all P + 10, incl cost
    176 0.002   CT_TAREALOAD    1   -CT_LOAD_DIS_PQ CT_ADD  20;     %% disp PQ + 10, incl cost
    177 0.002   CT_TAREALOAD    1   -CT_LOAD_DIS_P  CT_ADD  20;     %% disp P + 10, incl cost
];

%%-----  bus  -----
t = 'bus, replace';
g = apply_changes(1, mpc, chgtab);
e = mpc;
e.bus(1, GS) = 5;
t_is(g.bus, e.bus, 12, t);

t = 'bus, scale';
g = apply_changes(2, mpc, chgtab);
e = mpc;
e.bus(2, VMIN) = 1;
t_is(g.bus, e.bus, 12, t);

t = 'bus, shift';
g = apply_changes(3, mpc, chgtab);
e = mpc;
e.bus(:, VMAX) = mpc.bus(:, VMAX) - 0.01;
t_is(g.bus, e.bus, 12, t);

%%-----  branch  -----
t = 'branch, replace';
g = apply_changes(4, mpc, chgtab);
e = mpc;
e.branch(1, BR_STATUS) = 0;
t_is(g.branch, e.branch, 12, t);

t = 'branch, scale';
g = apply_changes(5, mpc, chgtab);
e = mpc;
e.branch(2, RATE_A) = 1000;
t_is(g.branch, e.branch, 12, t);

t = 'branch, shift';
g = apply_changes(6, mpc, chgtab);
e = mpc;
e.branch(:, RATE_B) = mpc.branch(:, RATE_B) + 0.1;
t_is(g.branch, e.branch, 12, t);

%%-----  gen  -----
t = 'gen, replace';
g = apply_changes(7, mpc, chgtab);
e = mpc;
e.gen(1, GEN_STATUS) = 0;
t_is(g.gen, e.gen, 12, t);

t = 'gen, scale';
g = apply_changes(8, mpc, chgtab);
e = mpc;
e.gen(2, QMAX) = 66;
t_is(g.gen, e.gen, 12, t);

t = 'gen, shift';
g = apply_changes(9, mpc, chgtab);
e = mpc;
e.gen(:, PMIN) = mpc.gen(:, PMIN) + 0.5;
t_is(g.gen, e.gen, 12, t);

%%-----  gencost  -----
t = 'gencost, replace';
g = apply_changes(78, mpc, chgtab);
e = mpc;
e.gencost(1, STARTUP) = 1000;
t_is(g.gencost, e.gencost, 12, t);

t = 'gencost, scale';
g = apply_changes(79, mpc, chgtab);
e = mpc;
e.gencost(2, COST+3) = 264;
t_is(g.gencost, e.gencost, 12, t);

t = 'gencost, shift';
g = apply_changes(80, mpc, chgtab);
e = mpc;
e.gencost(:, NCOST) = 3;
t_is(g.gencost, e.gencost, 12, t);

t = 'gencost, modcost x, scale';
g = apply_changes(81, mpc, chgtab);
e = mpc;
e.gencost(1, :) = modcost(mpc.gencost(1, :), 1.1, 'SCALE_X');
t_is(g.gencost, e.gencost, 12, t);

t = 'gencost, modcost x, scale all';
g = apply_changes(82, mpc, chgtab);
e = mpc;
e.gencost = modcost(mpc.gencost, 0.9, 'SCALE_X');
t_is(g.gencost, e.gencost, 12, t);

t = 'gencost, modcost x, shift';
g = apply_changes(83, mpc, chgtab);
e = mpc;
e.gencost(2, :) = modcost(mpc.gencost(2, :), -10, 'SHIFT_X');
t_is(g.gencost, e.gencost, 12, t);

t = 'gencost, modcost x, shift all';
g = apply_changes(84, mpc, chgtab);
e = mpc;
e.gencost = modcost(mpc.gencost, 20, 'SHIFT_X');
t_is(g.gencost, e.gencost, 12, t);

t = 'gencost, modcost f, scale';
g = apply_changes(85, mpc, chgtab);
e = mpc;
e.gencost(1, :) = modcost(mpc.gencost(1, :), 1.1, 'SCALE_F');
t_is(g.gencost, e.gencost, 12, t);

t = 'gencost, modcost f, scale all';
g = apply_changes(86, mpc, chgtab);
e = mpc;
e.gencost = modcost(mpc.gencost, 0.9, 'SCALE_F');
t_is(g.gencost, e.gencost, 12, t);

t = 'gencost, modcost f, shift';
g = apply_changes(87, mpc, chgtab);
e = mpc;
e.gencost(2, :) = modcost(mpc.gencost(2, :), -10, 'SHIFT_F');
t_is(g.gencost, e.gencost, 12, t);

t = 'gencost, modcost f, shift all';
g = apply_changes(88, mpc, chgtab);
e = mpc;
e.gencost = modcost(mpc.gencost, 20, 'SHIFT_F');
t_is(g.gencost, e.gencost, 12, t);

%%-----  area bus  -----
t = 'area bus : replace';
g = apply_changes(10, mpc, chgtab);
e = mpc;
e.bus(a{1}, VMAX) = 1.1;
t_is(g.bus, e.bus, 12, t);

t = 'area bus, scale';
g = apply_changes(11, mpc, chgtab);
e = mpc;
e.bus(a{2}, VMAX) = e.bus(a{2}, VMAX) * 1.01;
t_is(g.bus, e.bus, 12, t);

t = 'area bus, shift';
g = apply_changes(12, mpc, chgtab);
e = mpc;
e.bus(a{3}, VMAX) = mpc.bus(a{3}, VMAX) + 0.1;
t_is(g.bus, e.bus, 12, t);

%%-----  area branch  -----
t = 'area branch : replace';
g = apply_changes(13, mpc, chgtab);
e = mpc;
e.branch(bra{1}, RATE_B) = 100;
t_is(g.branch, e.branch, 12, t);

t = 'area branch, scale';
g = apply_changes(14, mpc, chgtab);
e = mpc;
e.branch(bra{2}, RATE_B) = e.branch(bra{2}, RATE_B) * 1.1;
t_is(g.branch, e.branch, 12, t);

t = 'area branch, shift';
g = apply_changes(15, mpc, chgtab);
e = mpc;
e.branch(bra{3}, RATE_B) = mpc.branch(bra{3}, RATE_B) + 0.5;
t_is(g.branch, e.branch, 12, t);

%%-----  area gen  -----
t = 'area gen : replace';
g = apply_changes(16, mpc, chgtab);
e = mpc;
e.gen(ga{1}, PMIN) = 0;
t_is(g.gen, e.gen, 12, t);

t = 'area gen, scale';
g = apply_changes(17, mpc, chgtab);
e = mpc;
e.gen(ga{2}, PMIN) = e.gen(ga{2}, PMIN) * 1.1;
t_is(g.gen, e.gen, 12, t);

t = 'area gen, shift';
g = apply_changes(18, mpc, chgtab);
e = mpc;
e.gen(ga{3}, PMIN) = mpc.gen(ga{3}, PMIN) + 0.1;
t_is(g.gen, e.gen, 12, t);

%%-----  area gencost  -----
t = 'area gencost : replace';
g = apply_changes(89, mpc, chgtab);
e = mpc;
e.gencost(ga{1}, SHUTDOWN) = 500;
t_is(g.gencost, e.gencost, 12, t);

t = 'area gencost, scale';
g = apply_changes(90, mpc, chgtab);
e = mpc;
e.gencost(ga{2}, COST+5) = e.gencost(ga{2}, COST+5) * 1.1;
t_is(g.gencost, e.gencost, 12, t);

t = 'area gencost, shift';
g = apply_changes(91, mpc, chgtab);
e = mpc;
e.gencost(ga{3}, NCOST) = 2;
t_is(g.gencost, e.gencost, 12, t);

t = 'area gencost, modcost x, scale';
g = apply_changes(92, mpc, chgtab);
e = mpc;
e.gencost(ga{1}, :) = modcost(mpc.gencost(ga{1}, :), 1.1, 'SCALE_X');
t_is(g.gencost, e.gencost, 12, t);

t = 'area gencost, modcost x, shift';
g = apply_changes(93, mpc, chgtab);
e = mpc;
e.gencost(ga{2}, :) = modcost(mpc.gencost(ga{2}, :), -10, 'SHIFT_X');
t_is(g.gencost, e.gencost, 12, t);

t = 'area gencost, modcost f, scale';
g = apply_changes(94, mpc, chgtab);
e = mpc;
e.gencost(ga{1}, :) = modcost(mpc.gencost(ga{1}, :), 1.1, 'SCALE_F');
t_is(g.gencost, e.gencost, 12, t);

t = 'area gencost, modcost f, shift';
g = apply_changes(95, mpc, chgtab);
e = mpc;
e.gencost(ga{2}, :) = modcost(mpc.gencost(ga{2}, :), -10, 'SHIFT_F');
t_is(g.gencost, e.gencost, 12, t);

%%-----  single load zone, one scale factor  -----
load = 2;
t = 'all fixed loads (PQ) * 2 : ';
e = mpc;
g = apply_changes(20, mpc, chgtab);
t_is(sum(g.bus(:, PD)), load*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), load*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all fixed loads (P) * 2 : ';
e = mpc;
g = apply_changes(21, mpc, chgtab);
t_is(sum(g.bus(:, PD)), load*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all loads (PQ) * 2 : ';
e = mpc;
g = apply_changes(22, mpc, chgtab);
t_is(sum(g.bus(:, PD)), load*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), load*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), load*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), load*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all loads (P) * 2 : ';
e = mpc;
g = apply_changes(23, mpc, chgtab);
t_is(sum(g.bus(:, PD)), load*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all disp loads (PQ) * 2 : ';
e = mpc;
g = apply_changes(24, mpc, chgtab);
t_is(sum(g.bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), load*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), load*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all disp loads (P) * 2 : ';
e = mpc;
g = apply_changes(25, mpc, chgtab);
t_is(sum(g.bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all loads+cost (PQ) * 2 : ';
e = mpc;
e.gencost(ld, :) = modcost(e.gencost(ld, :), load, 'SCALE_F');
e.gencost(ld, :) = modcost(e.gencost(ld, :), load, 'SCALE_X');
g = apply_changes(122, mpc, chgtab);
t_is(sum(g.bus(:, PD)), load*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), load*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), load*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), load*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all loads+cost (P) * 2 : ';
e = mpc;
e.gencost(ld, :) = modcost(e.gencost(ld, :), load, 'SCALE_F');
e.gencost(ld, :) = modcost(e.gencost(ld, :), load, 'SCALE_X');
g = apply_changes(123, mpc, chgtab);
t_is(sum(g.bus(:, PD)), load*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all disp loads+cost (PQ) * 2 : ';
e = mpc;
e.gencost(ld, :) = modcost(e.gencost(ld, :), load, 'SCALE_F');
e.gencost(ld, :) = modcost(e.gencost(ld, :), load, 'SCALE_X');
g = apply_changes(124, mpc, chgtab);
t_is(sum(g.bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), load*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), load*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all disp loads+cost (P) * 2 : ';
e = mpc;
e.gencost(ld, :) = modcost(e.gencost(ld, :), load, 'SCALE_F');
e.gencost(ld, :) = modcost(e.gencost(ld, :), load, 'SCALE_X');
g = apply_changes(125, mpc, chgtab);
t_is(sum(g.bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

%%-----  single load zone, one replace quantity  -----
load = 200;
t = 'all fixed loads (PQ) => total = 200 : ';
e = mpc;
g = apply_changes(26, mpc, chgtab);
t_is(sum(g.bus(:, PD)), load-total.disp.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), (load-total.disp.p)/total.fixed.p*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all fixed loads (P) => total = 200 : ';
e = mpc;
g = apply_changes(27, mpc, chgtab);
t_is(sum(g.bus(:, PD)), load-total.disp.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all loads (PQ) => total = 200 : ';
e = mpc;
g = apply_changes(28, mpc, chgtab);
t_is(sum(g.bus(:, PD)), load/total.both.p*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), load/total.both.p*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load/total.both.p*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), load/total.both.p*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), load/total.both.p*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all loads (P) => total = 200 : ';
e = mpc;
g = apply_changes(29, mpc, chgtab);
t_is(sum(g.bus(:, PD)), load/total.both.p*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load/total.both.p*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all disp loads (PQ) => total = 200 : ';
e = mpc;
g = apply_changes(30, mpc, chgtab);
t_is(sum(g.bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load-total.fixed.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), (load-total.fixed.p)/total.disp.p*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), (load-total.fixed.p)/total.disp.p*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all disp loads (P) => total = 200 : ';
e = mpc;
g = apply_changes(31, mpc, chgtab);
t_is(sum(g.bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load-total.fixed.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all loads+cost (PQ) => total = 200 : ';
e = mpc;
e.gencost(ld, :) = modcost(e.gencost(ld, :), load/total.both.p, 'SCALE_F');
e.gencost(ld, :) = modcost(e.gencost(ld, :), load/total.both.p, 'SCALE_X');
g = apply_changes(128, mpc, chgtab);
t_is(sum(g.bus(:, PD)), load/total.both.p*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), load/total.both.p*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load/total.both.p*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), load/total.both.p*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), load/total.both.p*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all loads+cost (P) => total = 200 : ';
e = mpc;
e.gencost(ld, :) = modcost(e.gencost(ld, :), load/total.both.p, 'SCALE_F');
e.gencost(ld, :) = modcost(e.gencost(ld, :), load/total.both.p, 'SCALE_X');
g = apply_changes(129, mpc, chgtab);
t_is(sum(g.bus(:, PD)), load/total.both.p*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load/total.both.p*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all disp loads+cost (PQ) => total = 200 : ';
e = mpc;
e.gencost(ld, :) = modcost(e.gencost(ld, :), (load-total.fixed.p)/total.disp.p, 'SCALE_F');
e.gencost(ld, :) = modcost(e.gencost(ld, :), (load-total.fixed.p)/total.disp.p, 'SCALE_X');
g = apply_changes(130, mpc, chgtab);
t_is(sum(g.bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load-total.fixed.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), (load-total.fixed.p)/total.disp.p*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), (load-total.fixed.p)/total.disp.p*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all disp loads+cost (P) => total = 200 : ';
e = mpc;
e.gencost(ld, :) = modcost(e.gencost(ld, :), (load-total.fixed.p)/total.disp.p, 'SCALE_F');
e.gencost(ld, :) = modcost(e.gencost(ld, :), (load-total.fixed.p)/total.disp.p, 'SCALE_X');
g = apply_changes(131, mpc, chgtab);
t_is(sum(g.bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load-total.fixed.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

%%-----  single load zone, one shift quantity  -----
load = 25 + total_load(mpc, 'all');
t = 'all fixed loads (PQ) => total + 25 : ';
e = mpc;
g = apply_changes(32, mpc, chgtab);
t_is(sum(g.bus(:, PD)), load-total.disp.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), (load-total.disp.p)/total.fixed.p*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all fixed loads (P) => total + 25 : ';
e = mpc;
g = apply_changes(33, mpc, chgtab);
t_is(sum(g.bus(:, PD)), load-total.disp.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all loads (PQ) => total + 25 : ';
e = mpc;
g = apply_changes(34, mpc, chgtab);
t_is(sum(g.bus(:, PD)), load/total.both.p*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), load/total.both.p*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load/total.both.p*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), load/total.both.p*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), load/total.both.p*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all loads (P) => total + 25 : ';
e = mpc;
g = apply_changes(35, mpc, chgtab);
t_is(sum(g.bus(:, PD)), load/total.both.p*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load/total.both.p*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all disp loads (PQ) => total + 25 : ';
e = mpc;
g = apply_changes(36, mpc, chgtab);
t_is(sum(g.bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load-total.fixed.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), (load-total.fixed.p)/total.disp.p*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), (load-total.fixed.p)/total.disp.p*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all disp loads (P) => total + 25 : ';
e = mpc;
g = apply_changes(37, mpc, chgtab);
t_is(sum(g.bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load-total.fixed.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all loads+cost (PQ) => total + 25 : ';
e = mpc;
e.gencost(ld, :) = modcost(e.gencost(ld, :), load/total.both.p, 'SCALE_F');
e.gencost(ld, :) = modcost(e.gencost(ld, :), load/total.both.p, 'SCALE_X');
g = apply_changes(134, mpc, chgtab);
t_is(sum(g.bus(:, PD)), load/total.both.p*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), load/total.both.p*total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load/total.both.p*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), load/total.both.p*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), load/total.both.p*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all loads+cost (P) => total + 25 : ';
e = mpc;
e.gencost(ld, :) = modcost(e.gencost(ld, :), load/total.both.p, 'SCALE_F');
e.gencost(ld, :) = modcost(e.gencost(ld, :), load/total.both.p, 'SCALE_X');
g = apply_changes(135, mpc, chgtab);
t_is(sum(g.bus(:, PD)), load/total.both.p*total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load/total.both.p*total.disp.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all disp loads+cost (PQ) => total + 25 : ';
e = mpc;
e.gencost(ld, :) = modcost(e.gencost(ld, :), (load-total.fixed.p)/total.disp.p, 'SCALE_F');
e.gencost(ld, :) = modcost(e.gencost(ld, :), (load-total.fixed.p)/total.disp.p, 'SCALE_X');
g = apply_changes(136, mpc, chgtab);
t_is(sum(g.bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load-total.fixed.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), (load-total.fixed.p)/total.disp.p*total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), (load-total.fixed.p)/total.disp.p*total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'all disp loads+cost (P) => total + 25 : ';
e = mpc;
e.gencost(ld, :) = modcost(e.gencost(ld, :), (load-total.fixed.p)/total.disp.p, 'SCALE_F');
e.gencost(ld, :) = modcost(e.gencost(ld, :), (load-total.fixed.p)/total.disp.p, 'SCALE_X');
g = apply_changes(137, mpc, chgtab);
t_is(sum(g.bus(:, PD)), total.fixed.p, 8, [t 'total fixed P']);
t_is(sum(g.bus(:, QD)), total.fixed.q, 8, [t 'total fixed Q']);
t_is(-sum(g.gen(ld, PMIN)), load-total.fixed.p, 8, [t 'total disp P']);
t_is(-sum(g.gen(ld, QMIN)), total.disp.qmin, 8, [t 'total disp Qmin']);
t_is(-sum(g.gen(ld, QMAX)), total.disp.qmax, 8, [t 'total disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

%%-----  single bus, one scale factor  -----
nb = size(mpc.bus, 1);
ng = size(mpc.gen, 1);
idx1 = [1 3:nb]';
idx2 = (2:ng)';
load = 2;
t = 'bus 2 fixed loads (PQ) * 2 : ';
e = mpc;
g = apply_changes(40, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), load*mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), load*mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 fixed loads (P) * 2 : ';
e = mpc;
g = apply_changes(41, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), load*mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 loads (PQ) * 2 : ';
e = mpc;
g = apply_changes(42, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), load*mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), load*mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), load*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), load*mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), load*mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 loads (P) * 2 : ';
e = mpc;
g = apply_changes(43, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), load*mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), load*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 disp loads (PQ) * 2 : ';
e = mpc;
g = apply_changes(44, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), load*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), load*mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), load*mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 disp loads (P) * 2 : ';
e = mpc;
g = apply_changes(45, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), load*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 loads+cost (PQ) * 2 : ';
e = mpc;
e.gencost(1, :) = modcost(e.gencost(1, :), load, 'SCALE_F');
e.gencost(1, :) = modcost(e.gencost(1, :), load, 'SCALE_X');
g = apply_changes(142, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), load*mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), load*mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), load*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), load*mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), load*mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 loads+cost (P) * 2 : ';
e = mpc;
e.gencost(1, :) = modcost(e.gencost(1, :), load, 'SCALE_F');
e.gencost(1, :) = modcost(e.gencost(1, :), load, 'SCALE_X');
g = apply_changes(143, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), load*mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), load*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 disp loads+cost (PQ) * 2 : ';
e = mpc;
e.gencost(1, :) = modcost(e.gencost(1, :), load, 'SCALE_F');
e.gencost(1, :) = modcost(e.gencost(1, :), load, 'SCALE_X');
g = apply_changes(144, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), load*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), load*mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), load*mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 disp loads+cost (P) * 2 : ';
e = mpc;
e.gencost(1, :) = modcost(e.gencost(1, :), load, 'SCALE_F');
e.gencost(1, :) = modcost(e.gencost(1, :), load, 'SCALE_X');
g = apply_changes(145, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), load*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

%%-----  single bus, one replace quantity  -----
load = 50;
t = 'bus 2 fixed loads (PQ) => Pd = 50 : ';
scale = (load+mpc.gen(1, PMIN)) / mpc.bus(2, PD);
e = mpc;
g = apply_changes(46, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), scale*mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), scale*mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 fixed loads (P) => Pd = 50 : ';
e = mpc;
g = apply_changes(47, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), scale*mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 loads (PQ) => Pd = 50 : ';
scale = load / (mpc.bus(2, PD) - mpc.gen(1, PMIN));
e = mpc;
g = apply_changes(48, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), scale*mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), scale*mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), scale*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), scale*mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), scale*mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 loads (P) => Pd = 50 : ';
e = mpc;
g = apply_changes(49, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), scale*mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), scale*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 disp loads (PQ) => Pd = 50 : ';
scale = (load - mpc.bus(2, PD)) / -mpc.gen(1, PMIN);
e = mpc;
g = apply_changes(50, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), scale*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), scale*mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), scale*mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 disp loads (P) => Pd = 50 : ';
e = mpc;
g = apply_changes(51, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), scale*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 loads+cost (PQ) => Pd = 50 : ';
scale = load / (mpc.bus(2, PD) - mpc.gen(1, PMIN));
e = mpc;
e.gencost(1, :) = modcost(e.gencost(1, :), scale, 'SCALE_F');
e.gencost(1, :) = modcost(e.gencost(1, :), scale, 'SCALE_X');
g = apply_changes(148, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), scale*mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), scale*mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), scale*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), scale*mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), scale*mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 loads+cost (P) => Pd = 50 : ';
e = mpc;
e.gencost(1, :) = modcost(e.gencost(1, :), scale, 'SCALE_F');
e.gencost(1, :) = modcost(e.gencost(1, :), scale, 'SCALE_X');
g = apply_changes(149, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), scale*mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), scale*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 disp loads+cost (PQ) => Pd = 50 : ';
scale = (load - mpc.bus(2, PD)) / -mpc.gen(1, PMIN);
e = mpc;
e.gencost(1, :) = modcost(e.gencost(1, :), scale, 'SCALE_F');
e.gencost(1, :) = modcost(e.gencost(1, :), scale, 'SCALE_X');
g = apply_changes(150, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), scale*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), scale*mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), scale*mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 disp loads+cost (P) => Pd = 50 : ';
e = mpc;
e.gencost(1, :) = modcost(e.gencost(1, :), scale, 'SCALE_F');
e.gencost(1, :) = modcost(e.gencost(1, :), scale, 'SCALE_X');
g = apply_changes(151, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), scale*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

%%-----  single bus, one shift quantity  -----
load = 10 + (mpc.bus(2, PD) - mpc.gen(1, PMIN));
t = 'bus 2 fixed loads (PQ) => total + 25 : ';
scale = (load+mpc.gen(1, PMIN)) / mpc.bus(2, PD);
e = mpc;
g = apply_changes(52, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), scale*mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), scale*mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 fixed loads (P) => total + 25 : ';
e = mpc;
g = apply_changes(53, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), scale*mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 loads (PQ) => total + 25 : ';
scale = load / (mpc.bus(2, PD) - mpc.gen(1, PMIN));
e = mpc;
g = apply_changes(54, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), scale*mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), scale*mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), scale*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), scale*mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), scale*mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 loads (P) => total + 25 : ';
e = mpc;
g = apply_changes(55, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), scale*mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), scale*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 disp loads (PQ) => total + 25 : ';
scale = (load - mpc.bus(2, PD)) / -mpc.gen(1, PMIN);
e = mpc;
g = apply_changes(56, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), scale*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), scale*mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), scale*mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 disp loads (P) => total + 25 : ';
e = mpc;
g = apply_changes(57, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), scale*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 loads+cost (PQ) => total + 25 : ';
scale = load / (mpc.bus(2, PD) - mpc.gen(1, PMIN));
e = mpc;
e.gencost(1, :) = modcost(e.gencost(1, :), scale, 'SCALE_F');
e.gencost(1, :) = modcost(e.gencost(1, :), scale, 'SCALE_X');
g = apply_changes(154, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), scale*mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), scale*mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), scale*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), scale*mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), scale*mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 loads+cost (P) => total + 25 : ';
e = mpc;
e.gencost(1, :) = modcost(e.gencost(1, :), scale, 'SCALE_F');
e.gencost(1, :) = modcost(e.gencost(1, :), scale, 'SCALE_X');
g = apply_changes(155, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), scale*mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), scale*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 disp loads+cost (PQ) => total + 25 : ';
scale = (load - mpc.bus(2, PD)) / -mpc.gen(1, PMIN);
e = mpc;
e.gencost(1, :) = modcost(e.gencost(1, :), scale, 'SCALE_F');
e.gencost(1, :) = modcost(e.gencost(1, :), scale, 'SCALE_X');
g = apply_changes(156, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), scale*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), scale*mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), scale*mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

t = 'bus 2 disp loads+cost (P) => total + 25 : ';
e = mpc;
e.gencost(1, :) = modcost(e.gencost(1, :), scale, 'SCALE_F');
e.gencost(1, :) = modcost(e.gencost(1, :), scale, 'SCALE_X');
g = apply_changes(157, mpc, chgtab);
t_is(g.bus(idx1, :), mpc.bus(idx1, :), 12, [t 'bus']);
t_is(g.gen(idx2, :), mpc.gen(idx2, :), 12, [t 'gen']);
t_is(g.bus(2, PD), mpc.bus(2, PD), 8, [t 'fixed P']);
t_is(g.bus(2, QD), mpc.bus(2, QD), 8, [t 'fixed Q']);
t_is(g.gen(1, PMIN), scale*mpc.gen(1, PMIN), 8, [t 'disp P']);
t_is(g.gen(1, QMIN), mpc.gen(1, QMIN), 8, [t 'disp Qmin']);
t_is(g.gen(1, QMAX), mpc.gen(1, QMAX), 8, [t 'disp Qmax']);
t_is(g.gencost(ld, :), e.gencost(ld, :), 8, [t 'disp gencost']);

%%-----  1 zone, area scale factor  -----
t = 'area fixed loads (PQ) * [3 1 1] : ';
load = [3 1 1];
e = mpc;
g = apply_changes(60, mpc, chgtab);
for k = 1:length(load)
    if k > 1
        t_is(g.bus(a{k}, :), mpc.bus(a{k}, :), 12, sprintf('%s area %d bus', t, k));
        t_is(g.gen(ga{k}, :), mpc.gen(ga{k}, :), 12, sprintf('%s area %d gen', t, k));
    else
        t_is(sum(g.bus(a{k}, PD)), load(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
        t_is(sum(g.bus(a{k}, QD)), load(k)*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
        t_is(-sum(g.gen(lda{k}, PMIN)), area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
        t_is(-sum(g.gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(g.gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

t = 'area fixed loads (P) * [3 1 1] : ';
e = mpc;
g = apply_changes(61, mpc, chgtab);
for k = 1:length(load)
    if k > 1
        t_is(g.bus(a{k}, :), mpc.bus(a{k}, :), 12, sprintf('%s area %d bus', t, k));
        t_is(g.gen(ga{k}, :), mpc.gen(ga{k}, :), 12, sprintf('%s area %d gen', t, k));
    else
        t_is(sum(g.bus(a{k}, PD)), load(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
        t_is(sum(g.bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
        t_is(-sum(g.gen(lda{k}, PMIN)), area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
        t_is(-sum(g.gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(g.gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

t = 'all area loads (PQ) * [3 1 1] : ';
e = mpc;
g = apply_changes(62, mpc, chgtab);
for k = 1:length(load)
    if k > 1
        t_is(g.bus(a{k}, :), mpc.bus(a{k}, :), 12, sprintf('%s area %d bus', t, k));
        t_is(g.gen(ga{k}, :), mpc.gen(ga{k}, :), 12, sprintf('%s area %d gen', t, k));
    else
        t_is(sum(g.bus(a{k}, PD)), load(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
        t_is(sum(g.bus(a{k}, QD)), load(k)*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
        t_is(-sum(g.gen(lda{k}, PMIN)), load(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
        t_is(-sum(g.gen(lda{k}, QMIN)), load(k)*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(g.gen(lda{k}, QMAX)), load(k)*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

t = 'all area loads (P) * [3 1 1] : ';
e = mpc;
g = apply_changes(63, mpc, chgtab);
for k = 1:length(load)
    if k > 1
        t_is(g.bus(a{k}, :), mpc.bus(a{k}, :), 12, sprintf('%s area %d bus', t, k));
        t_is(g.gen(ga{k}, :), mpc.gen(ga{k}, :), 12, sprintf('%s area %d gen', t, k));
    else
        t_is(sum(g.bus(a{k}, PD)), load(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
        t_is(sum(g.bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
        t_is(-sum(g.gen(lda{k}, PMIN)), load(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
        t_is(-sum(g.gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(g.gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

t = 'area disp loads (PQ) * [3 1 1] : ';
e = mpc;
g = apply_changes(64, mpc, chgtab);
for k = 1:length(load)
    if k > 1
        t_is(g.bus(a{k}, :), mpc.bus(a{k}, :), 12, sprintf('%s area %d bus', t, k));
        t_is(g.gen(ga{k}, :), mpc.gen(ga{k}, :), 12, sprintf('%s area %d gen', t, k));
    else
        t_is(sum(g.bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
        t_is(sum(g.bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
        t_is(-sum(g.gen(lda{k}, PMIN)), load(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
        t_is(-sum(g.gen(lda{k}, QMIN)), load(k)*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(g.gen(lda{k}, QMAX)), load(k)*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

t = 'area disp loads (P) * [3 1 1] : ';
e = mpc;
g = apply_changes(65, mpc, chgtab);
for k = 1:length(load)
    if k > 1
        t_is(g.bus(a{k}, :), mpc.bus(a{k}, :), 12, sprintf('%s area %d bus', t, k));
        t_is(g.gen(ga{k}, :), mpc.gen(ga{k}, :), 12, sprintf('%s area %d gen', t, k));
    else
        t_is(sum(g.bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
        t_is(sum(g.bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
        t_is(-sum(g.gen(lda{k}, PMIN)), load(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
        t_is(-sum(g.gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(g.gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

%%-----  1 zones, area target quantity  -----
t = 'area fixed loads (PQ) => total = 100 : ';
load = 100;
e = mpc;
g = apply_changes(66, mpc, chgtab);
scale = [(load-area(1).disp.p) / area(1).fixed.p 1 1];
for k = 1:length(load)
    if k > 1
        t_is(g.bus(a{k}, :), mpc.bus(a{k}, :), 12, sprintf('%s area %d bus', t, k));
        t_is(g.gen(ga{k}, :), mpc.gen(ga{k}, :), 12, sprintf('%s area %d gen', t, k));
    else
        t_is(sum(g.bus(a{k}, PD)), scale(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
        t_is(sum(g.bus(a{k}, QD)), scale(k)*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
        t_is(-sum(g.gen(lda{k}, PMIN)), area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
        t_is(-sum(g.gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(g.gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

t = 'area fixed loads (P) => total = 100 : ';
e = mpc;
g = apply_changes(67, mpc, chgtab);
for k = 1:length(load)
    if k > 1
        t_is(g.bus(a{k}, :), mpc.bus(a{k}, :), 12, sprintf('%s area %d bus', t, k));
        t_is(g.gen(ga{k}, :), mpc.gen(ga{k}, :), 12, sprintf('%s area %d gen', t, k));
    else
        t_is(sum(g.bus(a{k}, PD)), scale(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
        t_is(sum(g.bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
        t_is(-sum(g.gen(lda{k}, PMIN)), area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
        t_is(-sum(g.gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(g.gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

t = 'all area loads (PQ) => total = 100 : ';
e = mpc;
g = apply_changes(68, mpc, chgtab);
scale = [load / (area(1).fixed.p + area(1).disp.p) 1 1];
for k = 1:length(load)
    if k > 1
        t_is(g.bus(a{k}, :), mpc.bus(a{k}, :), 12, sprintf('%s area %d bus', t, k));
        t_is(g.gen(ga{k}, :), mpc.gen(ga{k}, :), 12, sprintf('%s area %d gen', t, k));
    else
        t_is(sum(g.bus(a{k}, PD)), scale(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
        t_is(sum(g.bus(a{k}, QD)), scale(k)*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
        t_is(-sum(g.gen(lda{k}, PMIN)), scale(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
        t_is(-sum(g.gen(lda{k}, QMIN)), scale(k)*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(g.gen(lda{k}, QMAX)), scale(k)*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

t = 'all area loads (P) => total = 100 : ';
e = mpc;
g = apply_changes(69, mpc, chgtab);
for k = 1:length(load)
    if k > 1
        t_is(g.bus(a{k}, :), mpc.bus(a{k}, :), 12, sprintf('%s area %d bus', t, k));
        t_is(g.gen(ga{k}, :), mpc.gen(ga{k}, :), 12, sprintf('%s area %d gen', t, k));
    else
        t_is(sum(g.bus(a{k}, PD)), scale(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
        t_is(sum(g.bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
        t_is(-sum(g.gen(lda{k}, PMIN)), scale(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
        t_is(-sum(g.gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(g.gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

t = 'area disp loads (PQ) => total = 100 : ';
e = mpc;
g = apply_changes(70, mpc, chgtab);
scale = [(load - area(1).fixed.p) / area(1).disp.p 1 1];
for k = 1:length(load)
    if k > 1
        t_is(g.bus(a{k}, :), mpc.bus(a{k}, :), 12, sprintf('%s area %d bus', t, k));
        t_is(g.gen(ga{k}, :), mpc.gen(ga{k}, :), 12, sprintf('%s area %d gen', t, k));
    else
        t_is(sum(g.bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
        t_is(sum(g.bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
        t_is(-sum(g.gen(lda{k}, PMIN)), scale(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
        t_is(-sum(g.gen(lda{k}, QMIN)), scale(k)*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(g.gen(lda{k}, QMAX)), scale(k)*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

t = 'area disp loads (P) => total = 100 : ';
e = mpc;
g = apply_changes(71, mpc, chgtab);
for k = 1:length(load)
    if k > 1
        t_is(g.bus(a{k}, :), mpc.bus(a{k}, :), 12, sprintf('%s area %d bus', t, k));
        t_is(g.gen(ga{k}, :), mpc.gen(ga{k}, :), 12, sprintf('%s area %d gen', t, k));
    else
        t_is(sum(g.bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
        t_is(sum(g.bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
        t_is(-sum(g.gen(lda{k}, PMIN)), scale(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
        t_is(-sum(g.gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(g.gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

%%-----  1 zones, area shift quantity  -----
t = 'area fixed loads (PQ) => total = total + 20 : ';
load = area(1).both.p + 20;
e = mpc;
g = apply_changes(72, mpc, chgtab);
scale = [(load-area(1).disp.p) / area(1).fixed.p 1 1];
for k = 1:length(load)
    if k > 1
        t_is(g.bus(a{k}, :), mpc.bus(a{k}, :), 12, sprintf('%s area %d bus', t, k));
        t_is(g.gen(ga{k}, :), mpc.gen(ga{k}, :), 12, sprintf('%s area %d gen', t, k));
    else
        t_is(sum(g.bus(a{k}, PD)), scale(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
        t_is(sum(g.bus(a{k}, QD)), scale(k)*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
        t_is(-sum(g.gen(lda{k}, PMIN)), area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
        t_is(-sum(g.gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(g.gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

t = 'area fixed loads (P) => total + 20 : ';
e = mpc;
g = apply_changes(73, mpc, chgtab);
for k = 1:length(load)
    if k > 1
        t_is(g.bus(a{k}, :), mpc.bus(a{k}, :), 12, sprintf('%s area %d bus', t, k));
        t_is(g.gen(ga{k}, :), mpc.gen(ga{k}, :), 12, sprintf('%s area %d gen', t, k));
    else
        t_is(sum(g.bus(a{k}, PD)), scale(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
        t_is(sum(g.bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
        t_is(-sum(g.gen(lda{k}, PMIN)), area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
        t_is(-sum(g.gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(g.gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

t = 'all area loads (PQ) => total + 20 : ';
e = mpc;
g = apply_changes(74, mpc, chgtab);
scale = [load / (area(1).fixed.p + area(1).disp.p) 1 1];
for k = 1:length(load)
    if k > 1
        t_is(g.bus(a{k}, :), mpc.bus(a{k}, :), 12, sprintf('%s area %d bus', t, k));
        t_is(g.gen(ga{k}, :), mpc.gen(ga{k}, :), 12, sprintf('%s area %d gen', t, k));
    else
        t_is(sum(g.bus(a{k}, PD)), scale(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
        t_is(sum(g.bus(a{k}, QD)), scale(k)*area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
        t_is(-sum(g.gen(lda{k}, PMIN)), scale(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
        t_is(-sum(g.gen(lda{k}, QMIN)), scale(k)*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(g.gen(lda{k}, QMAX)), scale(k)*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

t = 'all area loads (P) => total + 20 : ';
e = mpc;
g = apply_changes(75, mpc, chgtab);
for k = 1:length(load)
    if k > 1
        t_is(g.bus(a{k}, :), mpc.bus(a{k}, :), 12, sprintf('%s area %d bus', t, k));
        t_is(g.gen(ga{k}, :), mpc.gen(ga{k}, :), 12, sprintf('%s area %d gen', t, k));
    else
        t_is(sum(g.bus(a{k}, PD)), scale(k)*area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
        t_is(sum(g.bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
        t_is(-sum(g.gen(lda{k}, PMIN)), scale(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
        t_is(-sum(g.gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(g.gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

t = 'area disp loads (PQ) => total + 20 : ';
e = mpc;
g = apply_changes(76, mpc, chgtab);
scale = [(load - area(1).fixed.p) / area(1).disp.p 1 1];
for k = 1:length(load)
    if k > 1
        t_is(g.bus(a{k}, :), mpc.bus(a{k}, :), 12, sprintf('%s area %d bus', t, k));
        t_is(g.gen(ga{k}, :), mpc.gen(ga{k}, :), 12, sprintf('%s area %d gen', t, k));
    else
        t_is(sum(g.bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
        t_is(sum(g.bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
        t_is(-sum(g.gen(lda{k}, PMIN)), scale(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
        t_is(-sum(g.gen(lda{k}, QMIN)), scale(k)*area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(g.gen(lda{k}, QMAX)), scale(k)*area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

t = 'area disp loads (P) => total + 20 : ';
e = mpc;
g = apply_changes(77, mpc, chgtab);
for k = 1:length(load)
    if k > 1
        t_is(g.bus(a{k}, :), mpc.bus(a{k}, :), 12, sprintf('%s area %d bus', t, k));
        t_is(g.gen(ga{k}, :), mpc.gen(ga{k}, :), 12, sprintf('%s area %d gen', t, k));
    else
        t_is(sum(g.bus(a{k}, PD)), area(k).fixed.p, 8, sprintf('%s area %d fixed P', t, k));
        t_is(sum(g.bus(a{k}, QD)), area(k).fixed.q, 8, sprintf('%s area %d fixed Q', t, k));
        t_is(-sum(g.gen(lda{k}, PMIN)), scale(k)*area(k).disp.p, 8, sprintf('%s area %d disp P', t, k));
        t_is(-sum(g.gen(lda{k}, QMIN)), area(k).disp.qmin, 8, sprintf('%s area %d disp Qmin', t, k));
        t_is(-sum(g.gen(lda{k}, QMAX)), area(k).disp.qmax, 8, sprintf('%s area %d disp Qmax', t, k));
    end
end

t_end;
