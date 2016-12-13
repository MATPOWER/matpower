function t_hasPQcap(quiet)
%T_HASPQCAP  Tests for HASPQCAP.

%   MATPOWER
%   Copyright (c) 2005-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

t_begin(4, quiet);

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
gen = [
	1	10	0	10	-10	1	100	1	10	2	0	0	0	0	0	0	0	0	0	0	0;
	1	10	0	10	-10	1	100	1	10	2	0	20	0	12	0	2	0	0	0	0	0;
	1	10	0	10	-10	1	100	1	10	2	0	20	-15	12	-15	2	0	0	0	0	0;
	1	10	0	10	-10	1	100	1	10	2	0	20	-12	0	-2	0	0	0	0	0	0;
	1	10	0	10	-10	1	100	1	10	2	0	20	-12	15	-2	15	0	0	0	0	0;
	1	10	0	10	-10	1	100	1	10	2	0	20	-12	12	-2	2	0	0	0	0	0;
	1	10	0	10	-10	1	100	1	10	2	0	20	0	12	0	8	0	0	0	0	0;
	1	10	0	10	-10	1	100	1	10	2	0	20	-15	12	-15	8	0	0	0	0	0;
	1	10	0	10	-10	1	100	1	10	2	0	20	-12	0	-8	0	0	0	0	0	0;
	1	10	0	10	-10	1	100	1	10	2	0	20	-12	15	-8	15	0	0	0	0	0;
	1	10	0	10	-10	1	100	1	10	2	0	20	-12	12	-8	8	0	0	0	0	0;
	1	10	0	10	-10	1	100	1	10	2	0	20	0	0	-20	20	0	0	0	0	0;
	1	10	0	10	-10	1	100	1	10	2	0	20	0	0	-22	18	0	0	0	0	0;
	1	10	0	10	-10	1	100	1	10	2	0	20	0	0	-18	22	0	0	0	0	0;
	1	10	0	0	0	1	100	1	10	2	2	10	0	0	0	0	0	0	0	0	0;
];

t = 'hasPQcap(gen)';
t_is(hasPQcap(gen), [0;1;1;1;1;1;1;0;1;0;0;1;1;1;0], 12, t);

t = 'hasPQcap(gen, ''B'')';
t_is(hasPQcap(gen, 'B'), [0;1;1;1;1;1;1;0;1;0;0;1;1;1;0], 12, t);

t = 'hasPQcap(gen, ''U'')';
t_is(hasPQcap(gen, 'U'), [0;1;1;1;0;1;0;0;1;0;0;1;1;1;0], 12, t);

t = 'hasPQcap(gen, ''L'')';
t_is(hasPQcap(gen, 'L'), [0;1;0;1;1;1;1;0;0;0;0;1;1;1;0], 12, t);

t_end;
