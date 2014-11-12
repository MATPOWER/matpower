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
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

if nargin < 1
    quiet = 0;
end

t_begin(106, quiet);

%%-----  mpc = ext2int/int2ext(mpc)  -----
t = 'mpc = ext2int(mpc) : ';
mpce = loadcase('t_case_ext');
mpci = loadcase('t_case_int');
mpc = ext2int(mpce);
t_is(mpc.bus, mpci.bus, 12, [t 'bus']);
t_is(mpc.branch, mpci.branch, 12, [t 'branch']);
t_is(mpc.gen, mpci.gen, 12, [t 'gen']);
t_is(mpc.gencost, mpci.gencost, 12, [t 'gencost']);
t_is(mpc.A, mpci.A, 12, [t 'A']);
t_is(mpc.N, mpci.N, 12, [t 'N']);
t = 'mpc = ext2int(mpc) - repeat : ';
mpc = ext2int(mpc);
t_is(mpc.bus, mpci.bus, 12, [t 'bus']);
t_is(mpc.branch, mpci.branch, 12, [t 'branch']);
t_is(mpc.gen, mpci.gen, 12, [t 'gen']);
t_is(mpc.gencost, mpci.gencost, 12, [t 'gencost']);
t_is(mpc.A, mpci.A, 12, [t 'A']);
t_is(mpc.N, mpci.N, 12, [t 'N']);
t = 'mpc = int2ext(mpc) : ';
mpc = int2ext(mpc);
t_is(mpc.bus, mpce.bus, 12, [t 'bus']);
t_is(mpc.branch, mpce.branch, 12, [t 'branch']);
t_is(mpc.gen, mpce.gen, 12, [t 'gen']);
t_is(mpc.gencost, mpce.gencost, 12, [t 'gencost']);
t_is(mpc.A, mpce.A, 12, [t 'A']);
t_is(mpc.N, mpce.N, 12, [t 'N']);

%%-----  val = e2i_data/i2e_data(mpc, val, ...)  -----
t = 'val = e2i_data(mpc, val, ''bus'')';
mpc = ext2int(mpce);
got = e2i_data(mpc, mpce.xbus, 'bus');
ex = mpce.xbus;
ex(6, :) = [];
t_is(got, ex, 12, t);
t = 'val = i2e_data(mpc, val, oldval, ''bus'')';
tmp = ones(size(mpce.xbus));
tmp(6, :) = mpce.xbus(6, :);
got = i2e_data(mpc, ex, tmp, 'bus');
t_is(got, mpce.xbus, 12, t);

t = 'val = e2i_data(mpc, val, ''bus'', 2)';
got = e2i_data(mpc, mpce.xbus, 'bus', 2);
ex = mpce.xbus;
ex(:, 6) = [];
t_is(got, ex, 12, t);
t = 'val = i2e_data(mpc, val, oldval, ''bus'', 2)';
tmp = ones(size(mpce.xbus));
tmp(:, 6) = mpce.xbus(:, 6);
got = i2e_data(mpc, ex, tmp, 'bus', 2);
t_is(got, mpce.xbus, 12, t);

t = 'val = e2i_data(mpc, val, ''gen'')';
got = e2i_data(mpc, mpce.xgen, 'gen');
ex = mpce.xgen([4 2 1], :);
t_is(got, ex, 12, t);
t = 'val = i2e_data(mpc, val, oldval, ''gen'')';
tmp = ones(size(mpce.xgen));
tmp(3, :) = mpce.xgen(3, :);
got = i2e_data(mpc, ex, tmp, 'gen');
t_is(got, mpce.xgen, 12, t);

t = 'val = e2i_data(mpc, val, ''gen'', 2)';
got = e2i_data(mpc, mpce.xgen, 'gen', 2);
ex = mpce.xgen(:, [4 2 1]);
t_is(got, ex, 12, t);
t = 'val = i2e_data(mpc, val, oldval, ''gen'', 2)';
tmp = ones(size(mpce.xgen));
tmp(:, 3) = mpce.xgen(:, 3);
got = i2e_data(mpc, ex, tmp, 'gen', 2);
t_is(got, mpce.xgen, 12, t);

t = 'val = e2i_data(mpc, val, ''branch'')';
got = e2i_data(mpc, mpce.xbranch, 'branch');
ex = mpce.xbranch;
ex(7, :) = [];
t_is(got, ex, 12, t);
t = 'val = i2e_data(mpc, val, oldval, ''branch'')';
tmp = ones(size(mpce.xbranch));
tmp(7, :) = mpce.xbranch(7, :);
got = i2e_data(mpc, ex, tmp, 'branch');
t_is(got, mpce.xbranch, 12, t);

t = 'val = e2i_data(mpc, val, ''branch'', 2)';
got = e2i_data(mpc, mpce.xbranch, 'branch', 2);
ex = mpce.xbranch;
ex(:, 7) = [];
t_is(got, ex, 12, t);
t = 'val = i2e_data(mpc, val, oldval, ''branch'', 2)';
tmp = ones(size(mpce.xbranch));
tmp(:, 7) = mpce.xbranch(:, 7);
got = i2e_data(mpc, ex, tmp, 'branch', 2);
t_is(got, mpce.xbranch, 12, t);

t = 'val = e2i_data(mpc, val, {''branch'', ''gen'', ''bus''})';
got = e2i_data(mpc, mpce.xrows, {'branch', 'gen', 'bus'});
ex = [mpce.xbranch([1:6, 8:10], 1:4); mpce.xgen([4 2 1], :); mpce.xbus([1:5, 7:10], 1:4); -ones(2, 4)];
t_is(got, ex, 12, t);
t = 'val = i2e_data(mpc, val, oldval, {''branch'', ''gen'', ''bus''})';
tmp1 = ones(size(mpce.xbranch(:, 1:4)));
tmp1(7, 1:4) = mpce.xbranch(7, 1:4);
tmp2 = ones(size(mpce.xgen));
tmp2(3, :) = mpce.xgen(3, :);
tmp3 = ones(size(mpce.xbus(:, 1:4)));
tmp3(6, 1:4) = mpce.xbus(6, 1:4);
tmp = [tmp1; tmp2; tmp3];
got = i2e_data(mpc, ex, tmp, {'branch', 'gen', 'bus'});
t_is(got, mpce.xrows, 12, t);

t = 'val = e2i_data(mpc, val, {''branch'', ''gen'', ''bus''}, 2)';
got = e2i_data(mpc, mpce.xcols, {'branch', 'gen', 'bus'}, 2);
ex = [mpce.xbranch([1:6, 8:10], 1:4); mpce.xgen([4 2 1], :); mpce.xbus([1:5, 7:10], 1:4); -ones(2, 4)]';
t_is(got, ex, 12, t);
t = 'val = i2e_data(mpc, val, oldval, {''branch'', ''gen'', ''bus''}, 2)';
tmp1 = ones(size(mpce.xbranch(:, 1:4)));
tmp1(7, 1:4) = mpce.xbranch(7, 1:4);
tmp2 = ones(size(mpce.xgen));
tmp2(3, :) = mpce.xgen(3, :);
tmp3 = ones(size(mpce.xbus(:, 1:4)));
tmp3(6, 1:4) = mpce.xbus(6, 1:4);
tmp = [tmp1; tmp2; tmp3]';
got = i2e_data(mpc, ex, tmp, {'branch', 'gen', 'bus'}, 2);
t_is(got, mpce.xcols, 12, t);

%%-----  val = e2i_data/i2e_data(mpc, cell, ...)  -----
t = 'val = e2i_data(mpc, cell, ''bus'')';
mpc = ext2int(mpce);
got = e2i_data(mpc, mpce.strbus, 'bus');
ex = mpce.strbus;
ex(6, :) = [];
t_is(cellfun(@str2num, got), cellfun(@str2num, ex), 12, t);
t = 'val = i2e_data(mpc, cell, oldval, ''bus'')';
tmp = cell(size(mpce.strbus));
tmp(6, :) = mpce.strbus(6, :);
got = i2e_data(mpc, ex, tmp, 'bus');
t_is(cellfun(@str2num, got), cellfun(@str2num, mpce.strbus), 12, t);

t = 'val = e2i_data(mpc, cell, ''bus'', 2)';
got = e2i_data(mpc, mpce.strbus, 'bus', 2);
ex = mpce.strbus;
ex(:, 6) = [];
t_is(cellfun(@str2num, got), cellfun(@str2num, ex), 12, t);
t = 'val = i2e_data(mpc, cell, oldval, ''bus'', 2)';
tmp = cell(size(mpce.strbus));
tmp(:, 6) = mpce.strbus(:, 6);
got = i2e_data(mpc, ex, tmp, 'bus', 2);
t_is(cellfun(@str2num, got), cellfun(@str2num, mpce.strbus), 12, t);

t = 'val = e2i_data(mpc, cell, ''gen'')';
got = e2i_data(mpc, mpce.strgen, 'gen');
ex = mpce.strgen([4 2 1], :);
t_is(cellfun(@str2num, got), cellfun(@str2num, ex), 12, t);
t = 'val = i2e_data(mpc, cell, oldval, ''gen'')';
tmp = cell(size(mpce.strgen));
tmp(3, :) = mpce.strgen(3, :);
got = i2e_data(mpc, ex, tmp, 'gen');
t_is(cellfun(@str2num, got), cellfun(@str2num, mpce.strgen), 12, t);

t = 'val = e2i_data(mpc, cell, ''gen'', 2)';
got = e2i_data(mpc, mpce.strgen, 'gen', 2);
ex = mpce.strgen(:, [4 2 1]);
t_is(cellfun(@str2num, got), cellfun(@str2num, ex), 12, t);
t = 'val = i2e_data(mpc, cell, oldval, ''gen'', 2)';
tmp = cell(size(mpce.strgen));
tmp(:, 3) = mpce.strgen(:, 3);
got = i2e_data(mpc, ex, tmp, 'gen', 2);
t_is(cellfun(@str2num, got), cellfun(@str2num, mpce.strgen), 12, t);

t = 'val = e2i_data(mpc, cell, ''branch'')';
got = e2i_data(mpc, mpce.strbranch, 'branch');
ex = mpce.strbranch;
ex(7, :) = [];
t_is(cellfun(@str2num, got), cellfun(@str2num, ex), 12, t);
t = 'val = i2e_data(mpc, cell, oldval, ''branch'')';
tmp = cell(size(mpce.strbranch));
tmp(7, :) = mpce.strbranch(7, :);
got = i2e_data(mpc, ex, tmp, 'branch');
t_is(cellfun(@str2num, got), cellfun(@str2num, mpce.strbranch), 12, t);

t = 'val = e2i_data(mpc, cell, ''branch'', 2)';
got = e2i_data(mpc, mpce.strbranch, 'branch', 2);
ex = mpce.strbranch;
ex(:, 7) = [];
t_is(cellfun(@str2num, got), cellfun(@str2num, ex), 12, t);
t = 'val = i2e_data(mpc, cell, oldval, ''branch'', 2)';
tmp = cell(size(mpce.strbranch));
tmp(:, 7) = mpce.strbranch(:, 7);
got = i2e_data(mpc, ex, tmp, 'branch', 2);
t_is(cellfun(@str2num, got), cellfun(@str2num, mpce.strbranch), 12, t);

t = 'val = e2i_data(mpc, cell, {''branch'', ''gen'', ''bus''})';
got = e2i_data(mpc, mpce.strrows, {'branch', 'gen', 'bus'});
ex = [mpce.strbranch([1:6, 8:10], 1:4); mpce.strgen([4 2 1], :); mpce.strbus([1:5, 7:10], 1:4); cellfun(@num2str, num2cell(-ones(2, 4)), 'UniformOutput', 0)];
t_is(cellfun(@str2num, got), cellfun(@str2num, ex), 12, t);
t = 'val = i2e_data(mpc, cell, oldval, {''branch'', ''gen'', ''bus''})';
tmp1 = cell(size(mpce.strbranch(:, 1:4)));
tmp1(7, 1:4) = mpce.strbranch(7, 1:4);
tmp2 = cell(size(mpce.strgen));
tmp2(3, :) = mpce.strgen(3, :);
tmp3 = cell(size(mpce.strbus(:, 1:4)));
tmp3(6, 1:4) = mpce.strbus(6, 1:4);
tmp = [tmp1; tmp2; tmp3];
got = i2e_data(mpc, ex, tmp, {'branch', 'gen', 'bus'});
t_is(cellfun(@str2num, got), cellfun(@str2num, mpce.strrows), 12, t);

t = 'val = e2i_data(mpc, cell, {''branch'', ''gen'', ''bus''}, 2)';
got = e2i_data(mpc, mpce.strcols, {'branch', 'gen', 'bus'}, 2);
ex = [mpce.strbranch([1:6, 8:10], 1:4); mpce.strgen([4 2 1], :); mpce.strbus([1:5, 7:10], 1:4); cellfun(@num2str, num2cell(-ones(2, 4)), 'UniformOutput', 0)]';
t_is(cellfun(@str2num, got), cellfun(@str2num, ex), 12, t);
t = 'val = i2e_data(mpc, cell, oldval, {''branch'', ''gen'', ''bus''}, 2)';
tmp1 = cell(size(mpce.strbranch(:, 1:4)));
tmp1(7, 1:4) = mpce.strbranch(7, 1:4);
tmp2 = cell(size(mpce.strgen));
tmp2(3, :) = mpce.strgen(3, :);
tmp3 = cell(size(mpce.strbus(:, 1:4)));
tmp3(6, 1:4) = mpce.strbus(6, 1:4);
tmp = [tmp1; tmp2; tmp3]';
got = i2e_data(mpc, ex, tmp, {'branch', 'gen', 'bus'}, 2);
t_is(cellfun(@str2num, got), cellfun(@str2num, mpce.strcols), 12, t);

%%-----  mpc = e2i_field/i2e_field(mpc, field, ...)  -----
t = 'mpc = e2i_field(mpc, field, ''bus'')';
mpc = ext2int(mpce);
ex = mpce.xbus;
ex(6, :) = [];
got = e2i_field(mpc, 'xbus', 'bus');
t_is(got.xbus, ex, 12, t);
t = 'mpc = i2e_field(mpc, field, ''bus'')';
got = i2e_field(got, 'xbus', 'bus');
t_is(got.xbus, mpce.xbus, 12, t);

t = 'mpc = e2i_field(mpc, field, ''bus'', 2)';
ex = mpce.xbus;
ex(:, 6) = [];
got = e2i_field(mpc, 'xbus', 'bus', 2);
t_is(got.xbus, ex, 12, t);
t = 'mpc = i2e_field(mpc, field, ''bus'', 2)';
got = i2e_field(got, 'xbus', 'bus', 2);
t_is(got.xbus, mpce.xbus, 12, t);

t = 'mpc = e2i_field(mpc, field, ''gen'')';
ex = mpce.xgen([4 2 1], :);
got = e2i_field(mpc, 'xgen', 'gen');
t_is(got.xgen, ex, 12, t);
t = 'mpc = i2e_field(mpc, field, ''gen'')';
got = i2e_field(got, 'xgen', 'gen');
t_is(got.xgen, mpce.xgen, 12, t);

t = 'mpc = e2i_field(mpc, field, ''gen'', 2)';
ex = mpce.xgen(:, [4 2 1]);
got = e2i_field(mpc, 'xgen', 'gen', 2);
t_is(got.xgen, ex, 12, t);
t = 'mpc = i2e_field(mpc, field, ''gen'', 2)';
got = i2e_field(got, 'xgen', 'gen', 2);
t_is(got.xgen, mpce.xgen, 12, t);

t = 'mpc = e2i_field(mpc, field, ''branch'')';
ex = mpce.xbranch;
ex(7, :) = [];
got = e2i_field(mpc, 'xbranch', 'branch');
t_is(got.xbranch, ex, 12, t);
t = 'mpc = i2e_field(mpc, field, ''branch'')';
got = i2e_field(got, 'xbranch', 'branch');
t_is(got.xbranch, mpce.xbranch, 12, t);

t = 'mpc = e2i_field(mpc, field, ''branch'', 2)';
ex = mpce.xbranch;
ex(:, 7) = [];
got = e2i_field(mpc, 'xbranch', 'branch', 2);
t_is(got.xbranch, ex, 12, t);
t = 'mpc = i2e_field(mpc, field, ''branch'', 2)';
got = i2e_field(got, 'xbranch', 'branch', 2);
t_is(got.xbranch, mpce.xbranch, 12, t);

t = 'mpc = e2i_field(mpc, field, {''branch'', ''gen'', ''bus''})';
ex = [mpce.xbranch([1:6, 8:10], 1:4); mpce.xgen([4 2 1], :); mpce.xbus([1:5, 7:10], 1:4); -ones(2, 4)];
got = e2i_field(mpc, 'xrows', {'branch', 'gen', 'bus'});
t_is(got.xrows, ex, 12, t);
t = 'mpc = i2e_field(mpc, field, {''branch'', ''gen'', ''bus''})';
got = i2e_field(got, 'xrows', {'branch', 'gen', 'bus'});
t_is(got.xrows, mpce.xrows, 12, t);

t = 'mpc = e2i_field(mpc, field, {''branch'', ''gen'', ''bus''}, 2)';
ex = [mpce.xbranch([1:6, 8:10], 1:4); mpce.xgen([4 2 1], :); mpce.xbus([1:5, 7:10], 1:4); -ones(2, 4)]';
got = e2i_field(mpc, 'xcols', {'branch', 'gen', 'bus'}, 2);
t_is(got.xcols, ex, 12, t);
t = 'mpc = i2e_field(mpc, field, {''branch'', ''gen'', ''bus''})';
got = i2e_field(got, 'xcols', {'branch', 'gen', 'bus'}, 2);
t_is(got.xcols, mpce.xcols, 12, t);

t = 'mpc = e2i_field(mpc, {''field1'', ''field2''}, ordering)';
ex = mpce.x.more([4 2 1], :);
got = e2i_field(mpc, {'x', 'more'}, 'gen');
t_is(got.x.more, ex, 12, t);
t = 'mpc = i2e_field(mpc, {''field1'', ''field2''}, ordering)';
got = i2e_field(got, {'x', 'more'}, 'gen');
t_is(got.x.more, mpce.x.more, 12, t);

t = 'mpc = e2i_field(mpc, {''field1'', ''field2''}, ordering, 2)';
ex = mpce.x.more(:, [4 2 1]);
got = e2i_field(mpc, {'x', 'more'}, 'gen', 2);
t_is(got.x.more, ex, 12, t);
t = 'mpc = i2e_field(mpc, {''field1'', ''field2''}, ordering, 2)';
got = i2e_field(got, {'x', 'more'}, 'gen', 2);
t_is(got.x.more, mpce.x.more, 12, t);

%%-----  mpc = e2i_field/i2e_field(mpc, cellfield, ...)  -----
t = 'mpc = e2i_field(mpc, cellfield, ''bus'')';
mpc = ext2int(mpce);
ex = mpce.strbus;
ex(6, :) = [];
got = e2i_field(mpc, 'strbus', 'bus');
t_is(cellfun(@str2num, got.strbus), cellfun(@str2num, ex), 12, t);
t = 'mpc = i2e_field(mpc, cellfield, ''bus'')';
got = i2e_field(got, 'strbus', 'bus');
t_is(cellfun(@str2num, got.strbus), cellfun(@str2num, mpce.strbus), 12, t);

t = 'mpc = e2i_field(mpc, cellfield, ''bus'', 2)';
ex = mpce.strbus;
ex(:, 6) = [];
got = e2i_field(mpc, 'strbus', 'bus', 2);
t_is(cellfun(@str2num, got.strbus), cellfun(@str2num, ex), 12, t);
t = 'mpc = i2e_field(mpc, cellfield, ''bus'', 2)';
got = i2e_field(got, 'strbus', 'bus', 2);
t_is(cellfun(@str2num, got.strbus), cellfun(@str2num, mpce.strbus), 12, t);

t = 'mpc = e2i_field(mpc, cellfield, ''gen'')';
ex = mpce.strgen([4 2 1], :);
got = e2i_field(mpc, 'strgen', 'gen');
t_is(cellfun(@str2num, got.strgen), cellfun(@str2num, ex), 12, t);
t = 'mpc = i2e_field(mpc, cellfield, ''gen'')';
got = i2e_field(got, 'strgen', 'gen');
t_is(cellfun(@str2num, got.strgen), cellfun(@str2num, mpce.strgen), 12, t);

t = 'mpc = e2i_field(mpc, cellfield, ''gen'', 2)';
ex = mpce.strgen(:, [4 2 1]);
got = e2i_field(mpc, 'strgen', 'gen', 2);
t_is(cellfun(@str2num, got.strgen), cellfun(@str2num, ex), 12, t);
t = 'mpc = i2e_field(mpc, cellfield, ''gen'', 2)';
got = i2e_field(got, 'strgen', 'gen', 2);
t_is(cellfun(@str2num, got.strgen), cellfun(@str2num, mpce.strgen), 12, t);

t = 'mpc = e2i_field(mpc, cellfield, ''branch'')';
ex = mpce.strbranch;
ex(7, :) = [];
got = e2i_field(mpc, 'strbranch', 'branch');
t_is(cellfun(@str2num, got.strbranch), cellfun(@str2num, ex), 12, t);
t = 'mpc = i2e_field(mpc, cellfield, ''branch'')';
got = i2e_field(got, 'strbranch', 'branch');
t_is(cellfun(@str2num, got.strbranch), cellfun(@str2num, mpce.strbranch), 12, t);

t = 'mpc = e2i_field(mpc, cellfield, ''branch'', 2)';
ex = mpce.strbranch;
ex(:, 7) = [];
got = e2i_field(mpc, 'strbranch', 'branch', 2);
t_is(cellfun(@str2num, got.strbranch), cellfun(@str2num, ex), 12, t);
t = 'mpc = i2e_field(mpc, cellfield, ''branch'', 2)';
got = i2e_field(got, 'strbranch', 'branch', 2);
t_is(cellfun(@str2num, got.strbranch), cellfun(@str2num, mpce.strbranch), 12, t);

t = 'mpc = e2i_field(mpc, cellfield, {''branch'', ''gen'', ''bus''})';
ex = [mpce.strbranch([1:6, 8:10], 1:4); mpce.strgen([4 2 1], :); mpce.strbus([1:5, 7:10], 1:4); cellfun(@num2str, num2cell(-ones(2, 4)), 'UniformOutput', 0)];
got = e2i_field(mpc, 'strrows', {'branch', 'gen', 'bus'});
t_is(cellfun(@str2num, got.strrows), cellfun(@str2num, ex), 12, t);
t = 'mpc = i2e_field(mpc, cellfield, {''branch'', ''gen'', ''bus''})';
got = i2e_field(got, 'strrows', {'branch', 'gen', 'bus'});
t_is(cellfun(@str2num, got.strrows), cellfun(@str2num, mpce.strrows), 12, t);

t = 'mpc = e2i_field(mpc, cellfield, {''branch'', ''gen'', ''bus''}, 2)';
ex = [mpce.strbranch([1:6, 8:10], 1:4); mpce.strgen([4 2 1], :); mpce.strbus([1:5, 7:10], 1:4); cellfun(@num2str, num2cell(-ones(2, 4)), 'UniformOutput', 0)]';
got = e2i_field(mpc, 'strcols', {'branch', 'gen', 'bus'}, 2);
t_is(cellfun(@str2num, got.strcols), cellfun(@str2num, ex), 12, t);
t = 'mpc = i2e_field(mpc, cellfield, {''branch'', ''gen'', ''bus''})';
got = i2e_field(got, 'strcols', {'branch', 'gen', 'bus'}, 2);
t_is(cellfun(@str2num, got.strcols), cellfun(@str2num, mpce.strcols), 12, t);

%%-----  more mpc = ext2int/int2ext(mpc)  -----
t = 'mpc = ext2int(mpc) - bus/gen/branch only : ';
mpce = loadcase('t_case_ext');
mpci = loadcase('t_case_int');
mpce = rmfield(mpce, 'gencost');
mpce = rmfield(mpce, 'A');
mpce = rmfield(mpce, 'N');
mpci = rmfield(mpci, 'gencost');
mpci = rmfield(mpci, 'A');
mpci = rmfield(mpci, 'N');
mpc = ext2int(mpce);
t_is(mpc.bus, mpci.bus, 12, [t 'bus']);
t_is(mpc.branch, mpci.branch, 12, [t 'branch']);
t_is(mpc.gen, mpci.gen, 12, [t 'gen']);

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
t_is(mpc.A, mpci.A, 12, [t 'A']);
t_is(mpc.N, mpci.N, 12, [t 'N']);

t = 'mpc = int2ext(mpc) - A, N are DC sized : ';
mpc = int2ext(mpc);
t_is(mpc.bus, mpce.bus, 12, [t 'bus']);
t_is(mpc.branch, mpce.branch, 12, [t 'branch']);
t_is(mpc.gen, mpce.gen, 12, [t 'gen']);
t_is(mpc.gencost, mpce.gencost, 12, [t 'gencost']);
t_is(mpc.A, mpce.A, 12, [t 'A']);
t_is(mpc.N, mpce.N, 12, [t 'N']);

t_end;
