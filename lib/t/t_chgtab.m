function chgtab = t_chgtab
%T_CHGTAB  Returns a change table suitable for use with APPLY_CHANGES

define_constants;

%% Change Table
%	label	prob	table	row	col	chgtype	newval
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
