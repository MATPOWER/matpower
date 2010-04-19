function t_ext2int2ext(quiet)
%T_EXT2INT2EXT  Tests EXT2INT and INT2EXT.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009-2010 by Power System Engineering Research Center (PSERC)
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

t_begin(85, quiet);

if quiet
    verbose = 0;
else
    verbose = 1;
end

%%-----  mpc = ext2int/int2ext(mpc)  -----
t = 'mpc = ext2int(mpc) : ';
mpce = loadcase('t_case_ext');
mpci = loadcase('t_case_int');
mpc = ext2int(mpce);
t_is(mpc.bus, mpci.bus, 12, [t 'bus']);
t_is(mpc.branch, mpci.branch, 12, [t 'branch']);
t_is(mpc.gen, mpci.gen, 12, [t 'gen']);
t_is(mpc.gencost, mpci.gencost, 12, [t 'gencost']);
t_is(mpc.areas, mpci.areas, 12, [t 'areas']);
t_is(mpc.A, mpci.A, 12, [t 'A']);
t_is(mpc.N, mpci.N, 12, [t 'N']);
t = 'mpc = ext2int(mpc) - repeat : ';
mpc = ext2int(mpc);
t_is(mpc.bus, mpci.bus, 12, [t 'bus']);
t_is(mpc.branch, mpci.branch, 12, [t 'branch']);
t_is(mpc.gen, mpci.gen, 12, [t 'gen']);
t_is(mpc.gencost, mpci.gencost, 12, [t 'gencost']);
t_is(mpc.areas, mpci.areas, 12, [t 'areas']);
t_is(mpc.A, mpci.A, 12, [t 'A']);
t_is(mpc.N, mpci.N, 12, [t 'N']);
t = 'mpc = int2ext(mpc) : ';
mpc = int2ext(mpc);
t_is(mpc.bus, mpce.bus, 12, [t 'bus']);
t_is(mpc.branch, mpce.branch, 12, [t 'branch']);
t_is(mpc.gen, mpce.gen, 12, [t 'gen']);
t_is(mpc.gencost, mpce.gencost, 12, [t 'gencost']);
t_is(mpc.areas, mpce.areas, 12, [t 'areas']);
t_is(mpc.A, mpce.A, 12, [t 'A']);
t_is(mpc.N, mpce.N, 12, [t 'N']);

%%-----  val = ext2int/int2ext(mpc, val, ...)  -----
t = 'val = ext2int(mpc, val, ''bus'')';
mpc = ext2int(mpce);
got = ext2int(mpc, mpce.xbus, 'bus');
ex = mpce.xbus;
ex(6, :) = [];
t_is(got, ex, 12, t);
t = 'val = int2ext(mpc, val, oldval, ''bus'')';
tmp = ones(size(mpce.xbus));
tmp(6, :) = mpce.xbus(6, :);
got = int2ext(mpc, ex, tmp, 'bus');
t_is(got, mpce.xbus, 12, t);

t = 'val = ext2int(mpc, val, ''bus'', 2)';
got = ext2int(mpc, mpce.xbus, 'bus', 2);
ex = mpce.xbus;
ex(:, 6) = [];
t_is(got, ex, 12, t);
t = 'val = int2ext(mpc, val, oldval, ''bus'', 2)';
tmp = ones(size(mpce.xbus));
tmp(:, 6) = mpce.xbus(:, 6);
got = int2ext(mpc, ex, tmp, 'bus', 2);
t_is(got, mpce.xbus, 12, t);

t = 'val = ext2int(mpc, val, ''gen'')';
got = ext2int(mpc, mpce.xgen, 'gen');
ex = mpce.xgen([4 2 1], :);
t_is(got, ex, 12, t);
t = 'val = int2ext(mpc, val, oldval, ''gen'')';
tmp = ones(size(mpce.xgen));
tmp(3, :) = mpce.xgen(3, :);
got = int2ext(mpc, ex, tmp, 'gen');
t_is(got, mpce.xgen, 12, t);

t = 'val = ext2int(mpc, val, ''gen'', 2)';
got = ext2int(mpc, mpce.xgen, 'gen', 2);
ex = mpce.xgen(:, [4 2 1]);
t_is(got, ex, 12, t);
t = 'val = int2ext(mpc, val, oldval, ''gen'', 2)';
tmp = ones(size(mpce.xgen));
tmp(:, 3) = mpce.xgen(:, 3);
got = int2ext(mpc, ex, tmp, 'gen', 2);
t_is(got, mpce.xgen, 12, t);

t = 'val = ext2int(mpc, val, ''branch'')';
got = ext2int(mpc, mpce.xbranch, 'branch');
ex = mpce.xbranch;
ex(7, :) = [];
t_is(got, ex, 12, t);
t = 'val = int2ext(mpc, val, oldval, ''branch'')';
tmp = ones(size(mpce.xbranch));
tmp(7, :) = mpce.xbranch(7, :);
got = int2ext(mpc, ex, tmp, 'branch');
t_is(got, mpce.xbranch, 12, t);

t = 'val = ext2int(mpc, val, ''branch'', 2)';
got = ext2int(mpc, mpce.xbranch, 'branch', 2);
ex = mpce.xbranch;
ex(:, 7) = [];
t_is(got, ex, 12, t);
t = 'val = int2ext(mpc, val, oldval, ''branch'', 2)';
tmp = ones(size(mpce.xbranch));
tmp(:, 7) = mpce.xbranch(:, 7);
got = int2ext(mpc, ex, tmp, 'branch', 2);
t_is(got, mpce.xbranch, 12, t);

t = 'val = ext2int(mpc, val, {''branch'', ''gen'', ''bus''})';
got = ext2int(mpc, mpce.xrows, {'branch', 'gen', 'bus'});
ex = [mpce.xbranch([1:6, 8:10], 1:4); mpce.xgen([4 2 1], :); mpce.xbus([1:5, 7:10], 1:4); -ones(2, 4)];
t_is(got, ex, 12, t);
t = 'val = int2ext(mpc, val, oldval, {''branch'', ''gen'', ''bus''})';
tmp1 = ones(size(mpce.xbranch(:, 1:4)));
tmp1(7, 1:4) = mpce.xbranch(7, 1:4);
tmp2 = ones(size(mpce.xgen));
tmp2(3, :) = mpce.xgen(3, :);
tmp3 = ones(size(mpce.xbus(:, 1:4)));
tmp3(6, 1:4) = mpce.xbus(6, 1:4);
tmp = [tmp1; tmp2; tmp3];
got = int2ext(mpc, ex, tmp, {'branch', 'gen', 'bus'});
t_is(got, mpce.xrows, 12, t);

t = 'val = ext2int(mpc, val, {''branch'', ''gen'', ''bus''}, 2)';
got = ext2int(mpc, mpce.xcols, {'branch', 'gen', 'bus'}, 2);
ex = [mpce.xbranch([1:6, 8:10], 1:4); mpce.xgen([4 2 1], :); mpce.xbus([1:5, 7:10], 1:4); -ones(2, 4)]';
t_is(got, ex, 12, t);
t = 'val = int2ext(mpc, val, oldval, {''branch'', ''gen'', ''bus''}, 2)';
tmp1 = ones(size(mpce.xbranch(:, 1:4)));
tmp1(7, 1:4) = mpce.xbranch(7, 1:4);
tmp2 = ones(size(mpce.xgen));
tmp2(3, :) = mpce.xgen(3, :);
tmp3 = ones(size(mpce.xbus(:, 1:4)));
tmp3(6, 1:4) = mpce.xbus(6, 1:4);
tmp = [tmp1; tmp2; tmp3]';
got = int2ext(mpc, ex, tmp, {'branch', 'gen', 'bus'}, 2);
t_is(got, mpce.xcols, 12, t);

%%-----  mpc = ext2int/int2ext(mpc, field, ...)  -----
t = 'mpc = ext2int(mpc, field, ''bus'')';
mpc = ext2int(mpce);
ex = mpce.xbus;
ex(6, :) = [];
got = ext2int(mpc, 'xbus', 'bus');
t_is(got.xbus, ex, 12, t);
t = 'mpc = int2ext(mpc, field, ''bus'')';
got = int2ext(got, 'xbus', 'bus');
t_is(got.xbus, mpce.xbus, 12, t);

t = 'mpc = ext2int(mpc, field, ''bus'', 2)';
ex = mpce.xbus;
ex(:, 6) = [];
got = ext2int(mpc, 'xbus', 'bus', 2);
t_is(got.xbus, ex, 12, t);
t = 'mpc = int2ext(mpc, field, ''bus'', 2)';
got = int2ext(got, 'xbus', 'bus', 2);
t_is(got.xbus, mpce.xbus, 12, t);

t = 'mpc = ext2int(mpc, field, ''gen'')';
ex = mpce.xgen([4 2 1], :);
got = ext2int(mpc, 'xgen', 'gen');
t_is(got.xgen, ex, 12, t);
t = 'mpc = int2ext(mpc, field, ''gen'')';
got = int2ext(got, 'xgen', 'gen');
t_is(got.xgen, mpce.xgen, 12, t);

t = 'mpc = ext2int(mpc, field, ''gen'', 2)';
ex = mpce.xgen(:, [4 2 1]);
got = ext2int(mpc, 'xgen', 'gen', 2);
t_is(got.xgen, ex, 12, t);
t = 'mpc = int2ext(mpc, field, ''gen'', 2)';
got = int2ext(got, 'xgen', 'gen', 2);
t_is(got.xgen, mpce.xgen, 12, t);

t = 'mpc = ext2int(mpc, field, ''branch'')';
ex = mpce.xbranch;
ex(7, :) = [];
got = ext2int(mpc, 'xbranch', 'branch');
t_is(got.xbranch, ex, 12, t);
t = 'mpc = int2ext(mpc, field, ''branch'')';
got = int2ext(got, 'xbranch', 'branch');
t_is(got.xbranch, mpce.xbranch, 12, t);

t = 'mpc = ext2int(mpc, field, ''branch'', 2)';
ex = mpce.xbranch;
ex(:, 7) = [];
got = ext2int(mpc, 'xbranch', 'branch', 2);
t_is(got.xbranch, ex, 12, t);
t = 'mpc = int2ext(mpc, field, ''branch'', 2)';
got = int2ext(got, 'xbranch', 'branch', 2);
t_is(got.xbranch, mpce.xbranch, 12, t);

t = 'mpc = ext2int(mpc, field, {''branch'', ''gen'', ''bus''})';
ex = [mpce.xbranch([1:6, 8:10], 1:4); mpce.xgen([4 2 1], :); mpce.xbus([1:5, 7:10], 1:4); -ones(2, 4)];
got = ext2int(mpc, 'xrows', {'branch', 'gen', 'bus'});
t_is(got.xrows, ex, 12, t);
t = 'mpc = int2ext(mpc, field, {''branch'', ''gen'', ''bus''})';
got = int2ext(got, 'xrows', {'branch', 'gen', 'bus'});
t_is(got.xrows, mpce.xrows, 12, t);

t = 'mpc = ext2int(mpc, field, {''branch'', ''gen'', ''bus''}, 2)';
ex = [mpce.xbranch([1:6, 8:10], 1:4); mpce.xgen([4 2 1], :); mpce.xbus([1:5, 7:10], 1:4); -ones(2, 4)]';
got = ext2int(mpc, 'xcols', {'branch', 'gen', 'bus'}, 2);
t_is(got.xcols, ex, 12, t);
t = 'mpc = int2ext(mpc, field, {''branch'', ''gen'', ''bus''})';
got = int2ext(got, 'xcols', {'branch', 'gen', 'bus'}, 2);
t_is(got.xcols, mpce.xcols, 12, t);

t = 'mpc = ext2int(mpc, {''field1'', ''field2''}, ordering)';
ex = mpce.x.more([4 2 1], :);
got = ext2int(mpc, {'x', 'more'}, 'gen');
t_is(got.x.more, ex, 12, t);
t = 'mpc = int2ext(mpc, {''field1'', ''field2''}, ordering)';
got = int2ext(got, {'x', 'more'}, 'gen');
t_is(got.x.more, mpce.x.more, 12, t);

t = 'mpc = ext2int(mpc, {''field1'', ''field2''}, ordering, 2)';
ex = mpce.x.more(:, [4 2 1]);
got = ext2int(mpc, {'x', 'more'}, 'gen', 2);
t_is(got.x.more, ex, 12, t);
t = 'mpc = int2ext(mpc, {''field1'', ''field2''}, ordering, 2)';
got = int2ext(got, {'x', 'more'}, 'gen', 2);
t_is(got.x.more, mpce.x.more, 12, t);

%%-----  more mpc = ext2int/int2ext(mpc)  -----
t = 'mpc = ext2int(mpc) - bus/gen/branch only : ';
mpce = loadcase('t_case_ext');
mpci = loadcase('t_case_int');
mpce = rmfield(mpce, 'gencost');
mpce = rmfield(mpce, 'areas');
mpce = rmfield(mpce, 'A');
mpce = rmfield(mpce, 'N');
mpci = rmfield(mpci, 'gencost');
mpci = rmfield(mpci, 'areas');
mpci = rmfield(mpci, 'A');
mpci = rmfield(mpci, 'N');
mpc = ext2int(mpce);
t_is(mpc.bus, mpci.bus, 12, [t 'bus']);
t_is(mpc.branch, mpci.branch, 12, [t 'branch']);
t_is(mpc.gen, mpci.gen, 12, [t 'gen']);

t = 'mpc = ext2int(mpc) - no areas/A : ';
mpce = loadcase('t_case_ext');
mpci = loadcase('t_case_int');
mpce = rmfield(mpce, 'areas');
mpce = rmfield(mpce, 'A');
mpci = rmfield(mpci, 'areas');
mpci = rmfield(mpci, 'A');
mpc = ext2int(mpce);
t_is(mpc.bus, mpci.bus, 12, [t 'bus']);
t_is(mpc.branch, mpci.branch, 12, [t 'branch']);
t_is(mpc.gen, mpci.gen, 12, [t 'gen']);
t_is(mpc.gencost, mpci.gencost, 12, [t 'gencost']);
t_is(mpc.N, mpci.N, 12, [t 'N']);

t = 'mpc = ext2int(mpc) - Qg cost, no N : ';
mpce = loadcase('t_case_ext');
mpci = loadcase('t_case_int');
mpce = rmfield(mpce, 'N');
mpci = rmfield(mpci, 'N');
mpce.gencost = [mpce.gencost; mpce.gencost];
mpci.gencost = [mpci.gencost; mpci.gencost];
mpc = ext2int(mpce);
t_is(mpc.bus, mpci.bus, 12, [t 'bus']);
t_is(mpc.branch, mpci.branch, 12, [t 'branch']);
t_is(mpc.gen, mpci.gen, 12, [t 'gen']);
t_is(mpc.gencost, mpci.gencost, 12, [t 'gencost']);
t_is(mpc.areas, mpci.areas, 12, [t 'areas']);
t_is(mpc.A, mpci.A, 12, [t 'A']);

t = 'mpc = ext2int(mpc) - A, N are DC sized : ';
mpce = loadcase('t_case_ext');
mpci = loadcase('t_case_int');
eVmQgcols = [11:20 25:28]';
iVmQgcols = [10:18 22:24]';
mpce.A(:, eVmQgcols) = [];
mpce.N(:, eVmQgcols) = [];
mpci.A(:, iVmQgcols) = [];
mpci.N(:, iVmQgcols) = [];
mpc = ext2int(mpce);
t_is(mpc.bus, mpci.bus, 12, [t 'bus']);
t_is(mpc.branch, mpci.branch, 12, [t 'branch']);
t_is(mpc.gen, mpci.gen, 12, [t 'gen']);
t_is(mpc.gencost, mpci.gencost, 12, [t 'gencost']);
t_is(mpc.areas, mpci.areas, 12, [t 'areas']);
t_is(mpc.A, mpci.A, 12, [t 'A']);
t_is(mpc.N, mpci.N, 12, [t 'N']);
t = 'mpc = int2ext(mpc) - A, N are DC sized : ';
mpc = int2ext(mpc);
t_is(mpc.bus, mpce.bus, 12, [t 'bus']);
t_is(mpc.branch, mpce.branch, 12, [t 'branch']);
t_is(mpc.gen, mpce.gen, 12, [t 'gen']);
t_is(mpc.gencost, mpce.gencost, 12, [t 'gencost']);
t_is(mpc.areas, mpce.areas, 12, [t 'areas']);
t_is(mpc.A, mpce.A, 12, [t 'A']);
t_is(mpc.N, mpce.N, 12, [t 'N']);

t_end;
