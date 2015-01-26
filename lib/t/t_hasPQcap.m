function t_hasPQcap(quiet)
%T_HASPQCAP  Tests for HASPQCAP.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2005-2010 by Power System Engineering Research Center (PSERC)
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
];

t = 'hasPQcap(gen)';
t_is(hasPQcap(gen), [0;1;1;1;1;1;1;0;1;0;0;1;1;1], 12, t);

t = 'hasPQcap(gen, ''B'')';
t_is(hasPQcap(gen, 'B'), [0;1;1;1;1;1;1;0;1;0;0;1;1;1], 12, t);

t = 'hasPQcap(gen, ''U'')';
t_is(hasPQcap(gen, 'U'), [0;1;1;1;0;1;0;0;1;0;0;1;1;1], 12, t);

t = 'hasPQcap(gen, ''L'')';
t_is(hasPQcap(gen, 'L'), [0;1;0;1;1;1;1;0;0;0;0;1;1;1], 12, t);

t_end;
