function t_ext2int2ext(quiet)
%T_EXT2INT2EXT  Tests EXT2INT, INT2EXT, and related functions.
%   Includes tests for GET_REORDER, SET_REORDER, E2I_DATA, I2E_DATA
%   E2I_FIELD, I2E_FIELD, EXT2INT and INT2EXT.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

num_tests = 165;
t_begin(num_tests, quiet);

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

%% for MATLAB versions prior to R2012a (v 7.14)
if ~exist('isequaln')
    eval('isequaln = @isequalwithequalnans;');
end

if have_fcn('matlab', 'vnum') < 7.001
    t_skip(num_tests, 'test requires cellfun() construct not available before MATLAB 7.1');
else
    mpce = loadcase('t_case_ext');
    mpci = loadcase('t_case_int');

    An = mpce.xbus;
    As = mpce.strbus;

    %%-----  get_reorder  -----
    k = [3; 7; 4; 1];
    t = 'get_reorder(A, k, 1) : numeric';
    t_is(get_reorder(An, k, 1), An(k, :), 12, t);

    t = 'get_reorder(A, k, 2) : numeric';
    t_is(get_reorder(An, k, 2), An(:, k), 12, t);

    %%-----  set_reorder  -----
    k = (2:2:10)';
    t = 'set_reorder(A, B, k, 1) : numeric';
    B = An(k, :) * -1;
    got = set_reorder(An, B, k, 1);
    ex = An;
    ex(k, :) = -ex(k, :);
    t_is(got, ex, 12, t);

    t = 'set_reorder(A, B, k, 2) : numeric';
    B = An(:, k) * -1;
    got = set_reorder(An, B, k, 2);
    ex = An;
    ex(:, k) = -ex(:, k);
    t_is(got, ex, 12, t);

    t = 'set_reorder(A, Bshort, k, 1) : numeric';
    k = (2:2:8)';
    B = An(k, 1:8) * -1;
    got = set_reorder(An, B, k, 1);
    ex = An;
    ex(k, 1:8) = -ex(k, 1:8);
    t_is(got, ex, 12, t);

    t = 'set_reorder(A, Bshort, k, 2) : numeric';
    B = An(1:8, k) * -1;
    got = set_reorder(An, B, k, 2);
    ex = An;
    ex(1:8, k) = -ex(1:8, k);
    t_is(got, ex, 12, t);

    t = 'set_reorder(A, Blong, k, 1) : numeric';
    k = (2:2:10)';
    B = An(k, :) * -1;
    B = [B B];
    got = set_reorder(An, B, k, 1);
    ex = [An zeros(size(An))];
    ex(k, :) = [-An(k, :) -An(k, :)];
    t_is(got, ex, 12, t);

    t = 'set_reorder(A, Blong, k, 2) : numeric';
    B = An(:, k) * -1;
    B = [B; B];
    got = set_reorder(An, B, k, 2);
    ex = [An; zeros(size(An))];
    ex(:, k) = [-An(:, k); -An(:, k)];
    t_is(got, ex, 12, t);


    %%-----  get_reorder (cell)  -----
    k = [3; 7; 4; 1];
    t = 'get_reorder(A, k, 1) : cell';
    t_is(cellfun(@str2num, get_reorder(As, k, 1)), cellfun(@str2num, As(k, :)), 12, t);

    t = 'get_reorder(A, k, 2) : cell';
    t_is(cellfun(@str2num, get_reorder(As, k, 2)), cellfun(@str2num, As(:, k)), 12, t);

    %%-----  set_reorder (cell)  -----
    k = (2:2:10)';
    t = 'set_reorder(A, B, k, 1) : cell';
    B = cellfun(@num2str, num2cell(An(k, :) * -1), 'UniformOutput', 0);
    got = set_reorder(As, B, k, 1);
    ex = As;
    ex(k, :) = cellfun(@num2str, num2cell(-cellfun(@str2num, ex(k, :))), 'UniformOutput', 0);
    t_is(cellfun(@str2num, got), cellfun(@str2num, ex), 12, t);

    t = 'set_reorder(A, B, k, 2) : cell';
    B = cellfun(@num2str, num2cell(An(:, k) * -1), 'UniformOutput', 0);
    got = set_reorder(As, B, k, 2);
    ex = As;
    ex(:, k) = cellfun(@num2str, num2cell(-cellfun(@str2num, ex(:, k))), 'UniformOutput', 0);
    t_is(cellfun(@str2num, got), cellfun(@str2num, ex), 12, t);

    t = 'set_reorder(A, Bshort, k, 1) : cell';
    k = (2:2:8)';
    B = cellfun(@num2str, num2cell(An(k, 1:8) * -1), 'UniformOutput', 0);
    got = set_reorder(As, B, k, 1);
    ex = As;
    ex(k, 1:8) = cellfun(@num2str, num2cell(-cellfun(@str2num, ex(k, 1:8))), 'UniformOutput', 0);
    t_is(cellfun(@str2num, got), cellfun(@str2num, ex), 12, t);

    t = 'set_reorder(A, Bshort, k, 2) : cell';
    B = cellfun(@num2str, num2cell(An(1:8, k) * -1), 'UniformOutput', 0);
    got = set_reorder(As, B, k, 2);
    ex = As;
    ex(1:8, k) = cellfun(@num2str, num2cell(-cellfun(@str2num, ex(1:8, k))), 'UniformOutput', 0);
    t_is(cellfun(@str2num, got), cellfun(@str2num, ex), 12, t);

    t = 'set_reorder(A, Blong, k, 1) : cell';
    k = (2:2:10)';
    B = cellfun(@num2str, num2cell(An(k, :) * -1), 'UniformOutput', 0);
    B = [B B];
    got = set_reorder(As, B, k, 1);
    ex = [As cell(size(As))];
    ex(k, :) = cellfun(@num2str, num2cell([-An(k, :) -An(k, :)]), 'UniformOutput', 0);
    for i = 1:size(got, 1)      %% replace [] with '-999' to make str2num happy
        for j = 1:size(got, 2)
            if isempty(got{i, j})
                got{i, j} = '-999';
            end
            if isempty(ex{i, j})
                ex{i, j} = '-999';
            end
        end
    end
    t_is(cellfun(@str2num, got), cellfun(@str2num, ex), 12, t);

    t = 'set_reorder(A, Blong, k, 2) : cell';
    B = cellfun(@num2str, num2cell(An(:, k) * -1), 'UniformOutput', 0);
    B = [B; B];
    got = set_reorder(As, B, k, 2);
    ex = [As; cell(size(As))];
    ex(:, k) = cellfun(@num2str, num2cell([-An(:, k); -An(:, k)]), 'UniformOutput', 0);
    for i = 1:size(got, 1)      %% replace [] with '-999' to make str2num happy
        for j = 1:size(got, 2)
            if isempty(got{i, j})
                got{i, j} = '-999';
            end
            if isempty(ex{i, j})
                ex{i, j} = '-999';
            end
        end
    end
    t_is(cellfun(@str2num, got), cellfun(@str2num, ex), 12, t);

    %%-----  mpc = ext2int/int2ext(mpc)  -----
    t = 'mpc = ext2int(mpc) : ';
    mpc = ext2int(mpce);
    t_is(mpc.bus, mpci.bus, 12, [t 'bus']);
    t_is(mpc.branch, mpci.branch, 12, [t 'branch']);
    t_is(mpc.gen, mpci.gen, 12, [t 'gen']);
    t_is(mpc.gencost, mpci.gencost, 12, [t 'gencost']);
    t_ok(isequaln(mpc.bus_name, mpci.bus_name), [t 'bus_name']);
    t_ok(isequaln(mpc.gentype, mpci.gentype), [t 'gentype']);
    t_ok(isequaln(mpc.genfuel, mpci.genfuel), [t 'genfuel']);
    t_is(mpc.A, mpci.A, 12, [t 'A']);
    t_is(mpc.N, mpci.N, 12, [t 'N']);
    t = 'mpc = ext2int(mpc) - repeat : ';
    mpc = ext2int(mpc);
    t_is(mpc.bus, mpci.bus, 12, [t 'bus']);
    t_is(mpc.branch, mpci.branch, 12, [t 'branch']);
    t_is(mpc.gen, mpci.gen, 12, [t 'gen']);
    t_is(mpc.gencost, mpci.gencost, 12, [t 'gencost']);
    t_ok(isequaln(mpc.bus_name, mpci.bus_name), [t 'bus_name']);
    t_ok(isequaln(mpc.gentype, mpci.gentype), [t 'gentype']);
    t_ok(isequaln(mpc.genfuel, mpci.genfuel), [t 'genfuel']);
    t_is(mpc.A, mpci.A, 12, [t 'A']);
    t_is(mpc.N, mpci.N, 12, [t 'N']);
    t = 'mpc = int2ext(mpc) : ';
    mpc = int2ext(mpc);
    t_is(mpc.bus, mpce.bus, 12, [t 'bus']);
    t_is(mpc.branch, mpce.branch, 12, [t 'branch']);
    t_is(mpc.gen, mpce.gen, 12, [t 'gen']);
    t_is(mpc.gencost, mpce.gencost, 12, [t 'gencost']);
    t_ok(isequaln(mpc.bus_name, mpce.bus_name), [t 'bus_name']);
    t_ok(isequaln(mpc.gentype, mpce.gentype), [t 'gentype']);
    t_ok(isequaln(mpc.genfuel, mpce.genfuel), [t 'genfuel']);
    t_is(mpc.A, mpce.A, 12, [t 'A']);
    t_is(mpc.N, mpce.N, 12, [t 'N']);

    %%-----  mpc = ext2int/int2ext(mpc, mpopt)  -----
    mpopt = mpoption();
    t = 'mpc = ext2int(mpc, mpopt) : ';
    mpc = ext2int(mpce, mpopt);
    t_is(mpc.bus, mpci.bus, 12, [t 'bus']);
    t_is(mpc.branch, mpci.branch, 12, [t 'branch']);
    t_is(mpc.gen, mpci.gen, 12, [t 'gen']);
    t_is(mpc.gencost, mpci.gencost, 12, [t 'gencost']);
    t_ok(isequaln(mpc.bus_name, mpci.bus_name), [t 'bus_name']);
    t_ok(isequaln(mpc.gentype, mpci.gentype), [t 'gentype']);
    t_ok(isequaln(mpc.genfuel, mpci.genfuel), [t 'genfuel']);
    t_is(mpc.A, mpci.A, 12, [t 'A']);
    t_is(mpc.N, mpci.N, 12, [t 'N']);
    t = 'mpc = ext2int(mpc, mpopt) - repeat : ';
    mpc = ext2int(mpc, mpopt);
    t_is(mpc.bus, mpci.bus, 12, [t 'bus']);
    t_is(mpc.branch, mpci.branch, 12, [t 'branch']);
    t_is(mpc.gen, mpci.gen, 12, [t 'gen']);
    t_is(mpc.gencost, mpci.gencost, 12, [t 'gencost']);
    t_ok(isequaln(mpc.bus_name, mpci.bus_name), [t 'bus_name']);
    t_ok(isequaln(mpc.gentype, mpci.gentype), [t 'gentype']);
    t_ok(isequaln(mpc.genfuel, mpci.genfuel), [t 'genfuel']);
    t_is(mpc.A, mpci.A, 12, [t 'A']);
    t_is(mpc.N, mpci.N, 12, [t 'N']);
    t = 'mpc = int2ext(mpc, mpopt) : ';
    mpc = int2ext(mpc, mpopt);
    t_is(mpc.bus, mpce.bus, 12, [t 'bus']);
    t_is(mpc.branch, mpce.branch, 12, [t 'branch']);
    t_is(mpc.gen, mpce.gen, 12, [t 'gen']);
    t_is(mpc.gencost, mpce.gencost, 12, [t 'gencost']);
    t_ok(isequaln(mpc.bus_name, mpce.bus_name), [t 'bus_name']);
    t_ok(isequaln(mpc.gentype, mpce.gentype), [t 'gentype']);
    t_ok(isequaln(mpc.genfuel, mpce.genfuel), [t 'genfuel']);
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

    t = 'mpc = ext2int(mpc) - all buses isolated : ';
    mpc = loadcase('t_case_ext');
    mpc.bus(:, BUS_TYPE) = NONE;
    try
        mpci = ext2int(mpc);
        t_is(size(mpci.bus, 1), 0, 12, [t 'internal case empty']);
    catch
        t_ok(0, [t 'unexpected fatal error']);
    end

    t = 'mpc = int2ext(mpc) - all buses isolated : ';
    try
        mpce = int2ext(mpci);
        t_is(mpc.bus, mpce.bus, 12, [t 'bus']);
        t_is(mpc.branch, mpce.branch, 12, [t 'branch']);
        t_is(mpc.gen, mpce.gen, 12, [t 'gen']);
        t_is(mpc.gencost, mpce.gencost, 12, [t 'gencost']);
        t_is(mpc.A, mpce.A, 12, [t 'A']);
        t_is(mpc.N, mpce.N, 12, [t 'N']);
    catch
        t_ok(0, [t 'unexpected fatal error']);
        t_skip(5, [t 'unexpected fatal error']);
    end
end

t_end;
