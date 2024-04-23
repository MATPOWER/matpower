function obj = t_mp_table(quiet)
% t_mp_table - Tests for mp_table (and :class:`table`).

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

% define_constants;
if quiet
    verbose = 0;
else
    verbose = 1;
end

skip_tests_for_tablicious = 1;
table_classes = {@mp_table, @mp_table_subclass, @mp_foo_table};
class_names = {'mp_table', 'mp_table_subclass', 'mp_foo_table'};
if have_feature('table')
    table_classes = {@table, table_classes{:}};
    class_names = {'table', class_names{:}};
end
nc = length(table_classes);
nt = 229 + 8*nc;

t_begin(nc*nt, quiet);

%% set up example data
var1 = [1:6]';
var2 = var1 + var1/10;
var3 = {'one'; 'two'; 'three'; 'four'; 'five'; 'six'};
var4 = 1./var1;
var5 = [3;-5;0;1;-2;4] <= 0;
var6 = [10 15; 20 25; 30 35; 40 45; 50 55; 60 65];
var_names = {'igr', 'flt', 'str', 'dbl', 'boo', 'mat'};
var_values = {var1, var2, var3, var4, var5, var6};
row_names = {'row1'; 'second row'; 'row 3'; 'fourth row'; '5'; '6th row'};
dim_names = {'records', 'fields'};

%% run tests
for k = 1:nc
    table_class = table_classes{k};
    cls = upper(class_names{k});

    t = sprintf('%s : constructor - no args : ', cls);
    T = table_class();
    skip_oct_tab = isa(T, 'table') && have_feature('octave') && ...
            skip_tests_for_tablicious;
    skip_oct_tab2 = isa(T, 'mp_table_subclass') && have_feature('table') && ...
            have_feature('octave') && skip_tests_for_tablicious;
    skip_ml_tab = (isa(T, 'table') || isa(T, 'mp_table_subclass')) && ...
            have_feature('matlab') && have_feature('matlab', 'vnum') < 9.010;
    skip_ml_tab2 = (isa(T, 'table') || isa(T, 'mp_table_subclass')) && ...
            have_feature('matlab') && have_feature('matlab', 'vnum') < 9.012;
    t_ok(isa(T, class_names{k}), [t 'class']);
    t_ok(isempty(T), [t 'isempty']);

    t = sprintf('%s : constructor - ind vars : ', cls);
    T = table_class(var1, var2, var3, var4, var5(:), var6);
    t_ok(isa(T, class_names{k}), [t 'class']);
    t_ok(~isempty(T), [t 'not isempty']);
    t_is(size(T), [6 6], 12, [t 'sz = size(T)']);
    t_is(size(T, 1), 6, 12, [t 'sz = size(T, 1)']);
    t_is(size(T, 2), 6, 12, [t 'sz = size(T, 2)']);
    [nr, nz] = size(T);
    t_is([nr, nz], [6 6], 12, [t '[nr, nz] = size(T)']);
    if isa(T, 'mp_foo_table')
        t_skip(1, 'constructor does not capture var names');
    else
        t_ok(isequal(T.Properties.VariableNames, ...
            {'var1', 'var2', 'var3', 'var4', 'Var5', 'var6'}), [t 'VariableNames'] );
    end
    t_ok(isempty(T.Properties.RowNames), [t 'RowNames'] );
    if skip_ml_tab
        t_skip(1, [t 'DimensionNames not yet supported'] );
    else
        t_ok(isequal(T.Properties.DimensionNames, {'Row', 'Variables'}), ...
            [t 'DimensionNames'] );
    end
% show_me(T)

    %% subsref () (w/o RowNames)
    t = sprintf('%s : subsref () w/o RowNames : ', cls);
    t_is(size(T([2;4],[1;3;4])), [2 3], 12, [t 'T(ii,jj) : size(T(ii,jj))']);
    t_is(size(T(:,1:2:5)), [6 3], 12, [t 'T(:,j1:jN) : size']);
    t_is(size(T(1:2:5,:)), [3 6], 12, [t 'T(i1:iN,:) : size']);
    t_is(size(T(5,4)), [1 1], 12, [t 'T(i,j) : size']);
    t_is(size(T(4,6)), [1 1], 12, [t 'T(i,j) : size']);

    t = sprintf('%s : istable()', cls);
    t_ok(istable(T), t);

    t = sprintf('%s : isempty() : ', cls);
    t_ok(~isempty(T), [t 'false']);
    T1 = table_class();
    t_ok(isempty(T1), [t 'true']);

    t = sprintf('%s : constructor - ind vars, w/names : ', cls);
    if skip_ml_tab
        T = table_class(var1, var2, var3, var4, var5, var6, ...
            'VariableNames', var_names, 'RowNames', row_names);
    else
        T = table_class(var1, var2, var3, var4, var5, var6, ...
            'VariableNames', var_names, 'RowNames', row_names, ...
            'DimensionNames', dim_names);
    end
    t_ok(isa(T, class_names{k}), [t 'class'])
    t_ok(~isempty(T), [t 'not isempty']);
    t_is(size(T), [6 6], 12, [t 'sz = size(T)']);
    t_is(size(T, 1), 6, 12, [t 'sz = size(T, 1)']);
    t_is(size(T, 2), 6, 12, [t 'sz = size(T, 2)']);
    [nr, nz] = size(T);
    t_is([nr, nz], [6 6], 12, [t '[nr, nz] = size(T)']);
    t_ok(isequal(T.Properties.VariableNames, var_names), [t 'VariableNames'] );
    t_ok(isequal(T.Properties.RowNames, row_names), [t 'RowNames'] );
    if skip_ml_tab
        t_skip(1, [t 'DimensionNames not yet supported'] );
    else
        t_ok(isequal(T.Properties.DimensionNames, dim_names), [t 'DimensionNames'] );
    end

    %% . indexing
    %% get full variables
    t = sprintf('%s : subsref . : ', cls);
    t_is(T.igr, var1, 12, [t 'T.igr']);
    t_is(T.flt, var2, 12, [t 'T.flt']);
    t_ok(isequal(T.str, var3), [t 'T.str']);
    t_is(T.dbl, var4, 12, [t 'T.dbl']);
    t_is(T.boo, var5, 12, [t 'T.boo']);
    t_is(T.mat, var6, 12, [t 'T.mat']);

    %% get indexed variables
    t_is(T.igr(2), var1(2), 12, [t 'T.igr(i)']);
    t_is(T.flt(4:6), var2(4:6), 12, [t 'T.flt(i1:iN)']);
    t_ok(isequal(T.str{3}, var3{3}), [t 'T.str{i}']);
    t_ok(isequal(T.str(6:-1:4), var3(6:-1:4)), [t 'T.str(iN:-1:i1)']);
    t_is(T.dbl([5;3]), var4([5;3]), 12, [t 'T.dbl(ii)']);
    t_is(T.boo(var5 == 1), var5(var5 == 1), 12, [t 'T.boo(<logical>)']);
    t_is(T.mat([6;2;4], :), var6([6;2;4], :), 12, [t 'T.mat(ii, :)']);
    t_is(T.mat([6;2;4], 2), var6([6;2;4], 2), 12, [t 'T.mat(ii, j)']);
    t_is(T.mat([6;2;4], [2;1]), var6([6;2;4], [2;1]), 12, [t 'T.mat(ii, jj)']);
    t_is(T.mat(5, :), var6(5, :), 12, [t 'T.mat(i, :)']);
    t_is(T.mat(5, 2), var6(5, 2), 12, [t 'T.mat(i, j)']);
    t_is(T.mat(8:10), var6(8:10), 12, [t 'T.mat(jj)']);

    %% set full variables
    t = sprintf('%s : subsasgn . : ', cls);
    T.igr = var1(end:-1:1);
    T.flt = var2(end:-1:1);
    T.str = var3(end:-1:1);
    T.dbl = var4(end:-1:1);
    T.boo = var5(end:-1:1);
    T.mat = var6(end:-1:1, end:-1:1);
    t_is(T.igr, var1(end:-1:1), 12, [t 'T.igr']);
    t_is(T.flt, var2(end:-1:1), 12, [t 'T.flt']);
    t_ok(isequal(T.str, var3(end:-1:1)), [t 'T.str']);
    t_is(T.dbl, var4(end:-1:1), 12, [t 'T.dbl']);
    t_is(T.boo, var5(end:-1:1), 12, [t 'T.boo']);
    t_is(T.mat, var6(end:-1:1, end:-1:1), 12, [t 'T.mat']);

    %% set indexed variables
    T.igr(2) = 55;
    T.flt(4:6) = [0.4 0.5 0.6];
    T.str{3} = 'tres';
    T.str(6:-1:4) = {'seis'; 'cinco'; 'cuatro'};
    T.dbl([5;3]) = [pi; exp(1)];
    T.boo(var5 == 1) = false;
    T.mat([5;3],:) = [5 0.5; 3 0.3];

    v1 = var1(end:-1:1); v1(2) = 55;
    v2 = var2(end:-1:1); v2(4:6) = [0.4 0.5 0.6];
    v3 = var3(end:-1:1);
    v3(6:-1:3) = {'seis'; 'cinco'; 'cuatro'; 'tres'};
    v4 = var4(end:-1:1); v4([5;3]) = [pi; exp(1)];
    v5 = var5(end:-1:1); v5(var5 == 1) = false;
    v6 = var6(end:-1:1, end:-1:1); v6([5;3],:) = [5 0.5; 3 0.3];
    t_is(T.igr, v1, 12, [t 'T.igr(i)']);
    t_is(T.flt, v2, 12, [t 'T.flt(i1:iN)']);
    t_ok(isequal(T.str, v3), [t 'T.str{i}']);
    t_is(T.dbl, v4, 12, [t 'T.dbl(ii)']);
    t_is(T.boo, v5, 12, [t 'T.boo(<logical>)']);
    t_is(T.mat, v6, 12, [t 'T.mat(ii, :)']);

    %% {} indexing
    t = sprintf('%s : subsref {} : ', cls);
    t_ok(isequal(T{:, 1}, v1), [t 'T{:, 1} == v1']);
    t_ok(isequal(T{:, 2}, v2), [t 'T{:, 2} == v2']);
    t_ok(isequal(T{:, 3}, v3), [t 'T{:, 3} == v3']);
    t_ok(isequal(T{:, 4}, v4), [t 'T{:, 4} == v4']);
    t_ok(isequal(T{:, 5}, v5), [t 'T{:, 5} == v5']);
    t_ok(isequal(T{:, 6}, v6), [t 'T{:, 6} == v6']);

    t_ok(isequal(T{:, 4:6}, [v4 v5 v6]), [t 'T{:, j1:jN}']);
    t_ok(isequal(T{:, 2:2:6}, [v2 v4 v6]), [t 'T{:, j1:2:jN}']);

    t_ok(isequal(T{:, 'igr'}, v1), [t 'T{:, ''igr''} == v1']);
    if have_feature('matlab') && have_feature('matlab', 'vnum') < 9.005 % R2018b
        t_ok(isequal(T{:, 'flt'}, v2), [t 'T{:, ''flt''} == v2']);
    else
        cond = eval('isequal(T{:, "flt"}, v2)');
        t_ok(cond, [t 'T{:, "flt"} == v2']);
    end
    t_ok(isequal(T{:, 'str'}, v3), [t 'T{:, ''str''} == v3']);
    t_ok(isequal(T{:, 'dbl'}, v4), [t 'T{:, ''dbl''} == v4']);
    t_ok(isequal(T{:, 'boo'}, v5), [t 'T{:, ''boo''} == v5']);
    t_ok(isequal(T{:, 'mat'}, v6), [t 'T{:, ''mat''} == v6']);

    t_ok(isequal(T{2, 1}, v1(2)), [t 'T{i, j} == v<j>(i)']);
    t_ok(isequal(T{4, 6}, v6(4, :)), [t 'T{i, j} == v<j>(i)']);
    t_ok(isequal(T{end, end}, v6(end,:)), [t 'T{end, end} == v<end>(end,:)']);
    t_ok(isequal(T{1:3, 2}, v2(1:3)), [t 'T{i1:iN, j} == v<j>(i1:iN)']);
    t_ok(isequal(T{1:3, 6}, v6(1:3, :)), [t 'T{i1:iN, j} == v<j>(i1:iN, :)']);
    t_ok(isequal(T{[6;3], 5}, v5([6;3])), [t 'T{ii, j} == v<j>(ii)']);
    t_ok(isequal(T{[6;3], 6}, v6([6;3], :)), [t 'T{ii, j} == v<j>(ii, :)']);
    t_ok(isequal(T{6:-1:3, [2;4;5]}, [v2(6:-1:3) v4(6:-1:3) v5(6:-1:3)]), [t 'T{iN:-1:i1, jj}']);
    t_ok(isequal(T{6:-1:3, [2;4;6]}, [v2(6:-1:3) v4(6:-1:3) v6(6:-1:3, :)]), [t 'T{iN:-1:i1, jj}']);
    t_ok(isequal(T{:, [4;1]}, [v4 v1]), [t 'T{:, jj} == [v<j1> v<j2> ...]']);
    t_ok(isequal(T{:, [4;6;1]}, [v4 v6 v1]), [t 'T{:, jj} == [v<j1> v<j2> ...]']);

    t_ok(isequal(T{2, 'igr'}, v1(2)), [t 'T{i, <str>} == v<j>(i)']);
    t_ok(isequal(T{4, 'mat'}, v6(4, :)), [t 'T{i, <str>} == v<j>(i)']);
    t_ok(isequal(T{1:3, 'flt'}, v2(1:3)), [t 'T{i1:iN, <str>} == v<j>(i1:iN)']);
    t_ok(isequal(T{1:3, 'mat'}, v6(1:3, :)), [t 'T{i1:iN, <str>} == v<j>(i1:iN, :)']);
    t_ok(isequal(T{[6;3], 'boo'}, v5([6;3])), [t 'T{ii, <str>} == v<j>(ii)']);
    t_ok(isequal(T{[6;3], 'mat'}, v6([6;3], :)), [t 'T{ii, <str>} == v<j>(ii, :)']);
    t_ok(isequal(T{6:-1:3, {'flt', 'dbl', 'boo'}}, [v2(6:-1:3) v4(6:-1:3) v5(6:-1:3)]), [t 'T{iN:-1:i1, <cell>}']);
    t_ok(isequal(T{6:-1:3, {'flt', 'dbl', 'mat'}}, [v2(6:-1:3) v4(6:-1:3) v6(6:-1:3, :)]), [t 'T{iN:-1:i1, <cell>}']);
    t_ok(isequal(T{:, {'dbl', 'igr'}}, [v4 v1]), [t 'T{:, <cell>} == [v<j1> v<j2> ...]']);
    t_ok(isequal(T{:, {'dbl', 'mat', 'igr'}}, [v4 v6 v1]), [t 'T{:, <cell>} == [v<j1> v<j2> ...]']);

    t = sprintf('%s : subsasgn {} : ', cls);
    T{:, 1} = var1;
    t_is(T.igr, var1, 12, [t 'T{:, 1} = var1']);
    T{:, 2} = var2;
    t_is(T.flt, var2, 12, [t 'T{:, 2} = var2']);
    T{:, 3} = var3;
    t_ok(isequal(T.str, var3), [t 'T{:, 3} = var3']);
    T{:, 4} = var4;
    t_is(T.dbl, var4, 12, [t 'T{:, 4} = var4']);
    T{:, 5} = var5;
    t_is(T.boo, var5, 12, [t 'T{:, 5} = var5']);
    T{:, 6} = var6;
    t_is(T.mat, var6, 12, [t 'T{:, 6} = var6']);

    T{:, 4:6} = [v4 v5 v6];
    t_is(T{:, [2;4;5;6]}, [var2 v4 v5 v6], 12, [t 'T{:, j1:jN} = <mat>']);
    T{:, 2:2:6} = [v2 var4 var6];
    t_is(T{:, [2;4;5;6]}, [v2 var4 v5 var6], 12, [t 'T{:, j1:2:jN} = <mat>']);
    T{:, [2 5]} = [var2 var5];  %% reset values

    T2 = table_class([1;2;3], [4;5;6], [7;8;9]);
    T2{2, :} = T2{2, :} * 10;
    t_is(T2{:, :}, [1 4 7; 20 50 80; 3 6 9], 12, [t 'T{i, :} = <vec>']);
    T2{[1 3], :} = T2{[1 3], :} * -1;
    t_is(T2{:, :}, [-1 -4 -7; 20 50 80; -3 -6 -9], 12, [t 'T{ii, :} = <mat>']);
    T2{:, :} = [1 2 3; 4 5 6; 7 8 9];
    t_is(T2{:, :}, [1 2 3; 4 5 6; 7 8 9], 12, [t 'T{:, :} = <mat>']);

    T{:, 'igr'} = v1;
    t_is(T.igr, v1, 12, [t 'T{:, ''igr''} = v1']);
    T{:, 'flt'} = v2;
    t_is(T.flt, v2, 12, [t 'T{:, ''flt''} = v2']);
    T{:, 'str'} = v3;
    t_ok(isequal(T.str, v3), [t 'T{:, ''str''} = v3']);
    T{:, 'dbl'} = v4;
    t_is(T.dbl, v4, 12, [t 'T{:, ''dbl''} = v4']);
    T{:, 'boo'} = v5;
    t_is(T.boo, v5, 12, [t 'T{:, ''boo''} = v5']);
    T{:, 'mat'} = v6;
    t_is(T.mat, v6, 12, [t 'T{:, ''mat''} = v6']);
    T{:, {'str'}} = var3;
    t_ok(isequal(T{:, 3}, var3), [t 'T{:, <cell>} = var3']);
    T{:, {'igr', 'mat', 'dbl', 'flt', 'boo'}} = [var1 var6 var4 var2 var5];
    t_is(T{:, [1;2;4;5;6]}, [var1 var2 var4 var5 var6], 12, [t 'T{:, <cell>} = <mat>']);

    T{2, 1} = v1(2);
    t_ok(isequal(T{2, 1}, v1(2)), [t 'T{i, j} = v<j>(i)']);
    T{4, 6} = v6(4, :);
    t_ok(isequal(T{4, 6}, v6(4, :)), [t 'T{i, j} = v<j>(i, :)']);
    T{1:3, 2} = v2(1:3);
    t_ok(isequal(T.flt(1:3), v2(1:3)), [t 'T{i1:iN, j} = v<j>(i1:iN)']);
    T{1:3, 6} = v6(1:3, :);
    t_ok(isequal(T.mat(1:3, :), v6(1:3, :)), [t 'T{i1:iN, j} = v<j>(i1:iN, :)']);
    T{1:3, 6} = v6(1:3, [2; 1]);
    t_ok(isequal(T.mat(1:3, :), v6(1:3, [2; 1])), [t 'T{i1:iN, j} = v<j>(i1:iN, jj)']);
    T{[6;3], 5} = v5([6;3]);
    t_ok(isequal(T.boo([6;3]), v5([6;3])), [t 'T{ii, j} = v<j>(ii)']);
    T{[6;3], 6} = v6([6;3], :);
    t_ok(isequal(T.mat([6;3], :), v6([6;3], :)), [t 'T{ii, j} = v<j>(ii, :)']);
    T{[6;3], 6} = v6([6;3], [2;1]);
    t_ok(isequal(T.mat([6;3], :), v6([6;3], [2;1])), [t 'T{ii, j} = v<j>(ii, jj)']);
    T{6:-1:3, [2;4;6;5]} = [v2(6:-1:3) v4(6:-1:3) v6(6:-1:3, :) v5(6:-1:3)];
    t_ok(isequal(T{6:-1:3, [2;4;6;5]}, [v2(6:-1:3) v4(6:-1:3) v6(6:-1:3, :) v5(6:-1:3)]), [t 'T{iN:-1:i1, jj}']);
    T{:, [6;4;1]} = [v6 v4 v1];
    t_ok(isequal(T{:, [6;4;1]}, [v6 v4 v1]), [t 'T{:, jj} = [v<j1> v<j2> ...]']);

    T{2, 'igr'} = var1(2);
    t_ok(isequal(T{2, 1}, var1(2)), [t 'T{i, <str>} = var<j>(i)']);
    T{4, 'mat'} = var6(4, :);
    t_ok(isequal(T{4, 6}, var6(4, :)), [t 'T{i, <str>} = var<j>(i, :)']);
    T{1:3, 'flt'} = var2(1:3);
    t_ok(isequal(T.flt(1:3), var2(1:3)), [t 'T{i1:iN, <str>} = var<j>(i1:iN)']);
    T{1:3, 'mat'} = var6(1:3, :);
    t_ok(isequal(T.mat(1:3, :), var6(1:3, :)), [t 'T{i1:iN, <str>} = var<j>(i1:iN, :)']);
    T{1:3, 'mat'} = var6(1:3, [2; 1]);
    t_ok(isequal(T.mat(1:3, :), var6(1:3, [2; 1])), [t 'T{i1:iN, <str>} = var<j>(i1:iN, jj)']);
    T{[6;3], 'boo'} = var5([6;3]);
    t_ok(isequal(T.boo([6;3]), var5([6;3])), [t 'T{ii, <str>} = var<j>(ii)']);
    T{[6;3], 'mat'} = var6([6;3], :);
    t_ok(isequal(T.mat([6;3], :), var6([6;3], :)), [t 'T{ii, <str>} = var<j>(ii, :)']);
    T{[6;3], 'mat'} = var6([6;3], [2;1]);
    t_ok(isequal(T.mat([6;3], :), var6([6;3], [2;1])), [t 'T{ii, <str>} = v<j>(ii, jj)']);
    T{6:-1:3, {'flt', 'dbl', 'mat', 'boo'}} = [var2(6:-1:3) var4(6:-1:3) var6(6:-1:3, :) var5(6:-1:3)];
    t_ok(isequal(T{6:-1:3, [2;4;6;5]}, [var2(6:-1:3) var4(6:-1:3) var6(6:-1:3, :) var5(6:-1:3)]), [t 'T{iN:-1:i1, <cell>}']);
    T{:, {'mat', 'dbl', 'igr'}} = [var6 var4 var1];
    t_ok(isequal(T{:, [6;4;1]}, [var6 var4 var1]), [t 'T{:, <cell>} = [var<j1> var<j2> ...]']);

    T{:, 1} = v1;
    T{:, 2} = v2;
    T{:, 3} = v3;
    T{:, 4} = v4;
    T{:, 5} = v5;
    T{:, 6} = v6;

    %% () indexing
    t = sprintf('%s : subsref () : T(ii,jj) : ', cls);
    ii = [2;4]; jj = [1;3;4;6];
    T2 = T(ii,jj);
    t_is(size(T2), [2 4], 12, [t 'size']);
    t_ok(isequal(T2.Properties.VariableNames, var_names(jj)), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {v1(ii), v3(ii), v4(ii), v6(ii, :)}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names(ii)), [t 'RowNames']);
    if skip_ml_tab
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end

    t = sprintf('%s : subsref () : T(ii,<cell>) : ', cls);
    ii = [2;4]; jj = [1;3;4;6]; cc = {'igr', 'str', 'dbl', 'mat'};
    T2 = T(ii,cc);
    t_is(size(T2), [2 4], 12, [t 'size']);
    t_ok(isequal(T2.Properties.VariableNames, var_names(jj)), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {v1(ii), v3(ii), v4(ii), v6(ii, :)}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names(ii)), [t 'RowNames']);
    if skip_ml_tab
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end

    t = sprintf('%s : subsref () : T(:,j1:s:jN) : ', cls);
    jj = [2:2:6];
    T2 = T(:,2:2:6);
    t_is(size(T2), [6 3], 12, [t 'size']);
    t_ok(isequal(T2.Properties.VariableNames, var_names(jj)), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {v2, v4, v6}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names), [t 'RowNames']);
    if skip_ml_tab
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end

    t = sprintf('%s : subsref () : T(i1:s:iN,:)) : ', cls);
    ii = [1:2:5];
    T2 = T(1:2:5,:);
    t_is(size(T2), [3 6], 12, [t 'size']);
    t_ok(isequal(T2.Properties.VariableNames, var_names), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {v1(ii), v2(ii), v3(ii), v4(ii), v5(ii), v6(ii, :)}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names(ii)), [t 'RowNames']);
    if skip_ml_tab
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end

    t = sprintf('%s : subsref () : T(i,j) : ', cls);
    T2 = T(5,4);
    t_is(size(T2), [1 1], 12, [t 'size']);
    t_ok(isequal(T2.Properties.VariableNames, var_names(4)), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {v4(5)}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names(5)), [t 'RowNames']);
    if skip_ml_tab
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end
    T2 = T(3,6);
    t_is(size(T2), [1 1], 12, [t 'size']);
    t_ok(isequal(T2.Properties.VariableNames, var_names(6)), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {v6(3, :)}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names(3)), [t 'RowNames']);
    if skip_ml_tab
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end

    t = sprintf('%s : subsref () : T(i,<str>) : ', cls);
    T2 = T(5,'dbl');
    t_is(size(T2), [1 1], 12, [t 'size']);
    t_ok(isequal(T2.Properties.VariableNames, var_names(4)), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {v4(5)}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names(5)), [t 'RowNames']);
    if skip_ml_tab
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end
    T2 = T(3,'mat');
    t_is(size(T2), [1 1], 12, [t 'size']);
    t_ok(isequal(T2.Properties.VariableNames, var_names(6)), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {v6(3, :)}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names(3)), [t 'RowNames']);
    if skip_ml_tab
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end

    t = sprintf('%s : subsref () : T(end, end) : ', cls);
    T2 = T(end, end);
    t_is(size(T2), [1 1], 12, [t 'size']);
    t_ok(isequal(T2.Properties.VariableNames, var_names(6)), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {v6(end, :)}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names(end)), [t 'RowNames']);
    if skip_ml_tab
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end

    t = sprintf('%s : subsasgn () : T(ii,jj) : ', cls);
    ii0 = [2;4]; jj = [1;3;4;6];
    T3 = T(ii0,jj);
    ii = [4;2];
    T(ii,jj) = T3;
    T2 = T(ii,jj);
    t_ok(isequal(T2.Properties.VariableNames, var_names(jj)), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {v1(ii0), v3(ii0), v4(ii0), v6(ii0, :)}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names(ii)), [t 'RowNames']);
    if skip_ml_tab
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end
    T(ii0, jj) = T2;        %% restore

    t = sprintf('%s : subsasgn () : T(ii,<cell>) : ', cls);
    ii0 = [2;4]; jj = [1;3;4;6]; cc = {'igr', 'str', 'dbl', 'mat'};
    T3 = T(ii0,jj);
    ii = [4;2];
    T(ii,cc) = T3;
    T2 = T(ii,jj);
    t_ok(isequal(T2.Properties.VariableNames, var_names(jj)), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {v1(ii0), v3(ii0), v4(ii0), v6(ii0, :)}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names(ii)), [t 'RowNames']);
    if skip_ml_tab
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end
    T(ii0, jj) = T2;        %% restore

    t = sprintf('%s : subsasgn () : T(:,j1:s:jN) : ', cls);
    jj = [2:2:6];
    T3 = T(end:-1:1,2:2:6);
    T(:,2:2:6) = T3;
    T2 = T(:,2:2:6);
    t_ok(isequal(T2.Properties.VariableNames, var_names(jj)), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {v2(end:-1:1), v4(end:-1:1), v6(end:-1:1, :)}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names), [t 'RowNames']);
    if skip_ml_tab
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end
    T(end:-1:1,2:2:6) = T3;     %% restore

    t = sprintf('%s : subsasgn () : T(i1:s:iN,:)) : ', cls);
    ii0 = [1:2:5];
    T3 = T(1:2:5,:);
    ii = [5:-2:1];
    T(5:-2:1, :) = T3;
    T2 = T(5:-2:1, :);
    t_ok(isequal(T2.Properties.VariableNames, var_names), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {v1(ii0), v2(ii0), v3(ii0), v4(ii0), v5(ii0), v6(ii0, :)}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names(ii)), [t 'RowNames']);
    if skip_ml_tab
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end
    T(1:2:5, :) = T3;   %% restore

    t = sprintf('%s : subsref () : T(i,j) : ', cls);
    T3 = T(5,4);
    T3{1,1} = exp(1);
    T(5,4) = T3;
    T2 = T(5,4);
    t_ok(isequal(T2.Properties.VariableNames, var_names(4)), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {exp(1)}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names(5)), [t 'RowNames']);
    if skip_ml_tab
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end
    T{5,4} = pi;

    T3 = T(3,6);
    T3{1,1} = [111 222];
    T(3,6) = T3;
    T2 = T(3,6);
    t_ok(isequal(T2.Properties.VariableNames, var_names(6)), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {[111 222]}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names(3)), [t 'RowNames']);
    if skip_ml_tab
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end
    T{3,6} = [3 0.3];

    t = sprintf('%s : subsref () : T(i,<str>) : ', cls);
    T3 = T(5,4);
    T3{1,1} = exp(1);
    T(5,4) = T3;
    T2 = T(5,'dbl');
    t_ok(isequal(T2.Properties.VariableNames, var_names(4)), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {exp(1)}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names(5)), [t 'RowNames']);
    if skip_ml_tab
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end
    T{5,4} = pi;

    T3 = T(3,6);
    T3{1,1} = [111 222];
    T(3,6) = T3;
    T2 = T(3,'mat');
    t_ok(isequal(T2.Properties.VariableNames, var_names(6)), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {[111 222]}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names(3)), [t 'RowNames']);
    if skip_ml_tab
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end
    T{3,6} = [3 0.3];

    %% vertical concatenation
    t = sprintf('%s : vertical concatenation : ', cls);
    T1 = T(1:3, :);
    T2 = T(4:6, :);
    T3 = T([4:6 1:3], :);
    T4 = [T2; T1];
    t_ok(isequal(T4, T3), [t '[T1;T2]']);
    T5 = table_class([7;8],[7.7;8.8],{'seven';'eight'},1./[7;8],[-1;2]<=0, [70 75; 80 85], ...
        'VariableNames', var_names);
    if skip_oct_tab || skip_oct_tab2
        t_skip(1, [t 'RowNames auto-generation not yet supported']);
    else
        T6 = [T5; T2; T1];
        t_ok(isequal(T6, [T5;T3]), [t '[T1;T2;T3]']);
    end

    %% horizontal concatenation
    t = sprintf('%s : horizontal concatenation : ', cls);
    T1 = T(:, 1:3);
    T2 = T(:, 4:6);
    T3 = T(:, [4:6 1:3]);
    T4 = [T2 T1];
    t_ok(isequal(T4, T3), [t '[T1 T2]']);
    T5 = table_class(var3, var2);
    T6 = [T2 T1 T5];
    t_ok(isequal(T6, [T4 T5]), [t '[T1 T2 T3]']);

    %% deleting variables
    t = sprintf('%s : delete variables : ', cls);
    if skip_ml_tab2
        t_skip(1, [t 'not yet supported']);
    else
        T7 = T;
        T7(:, [5 6 4]) = [];
        t_ok(isequal(T7, T1), [t 'T(:, j) = []']);
    end

    %% adding variables
    t = sprintf('%s : add variables : ', cls);
    T7 = T1;
    T7.dbl = T2.dbl;
    T7.boo = T2.boo;
    T7.mat = T2.mat;
    t_ok(isequal(T7, T), [t 'T.new_var = val']);

    %% more {} indexing
    t = sprintf('%s : more subsref {} : ', cls);
    T2 = T(:, [1;2;4;6;5]);
    T3 = T6(:, 6:7);
    ii = [5;3;1];
    jj = [6:7];
    ex = horzcat(v3, var3);
    t_is(T2{ii, :}, [v1(ii) v2(ii) v4(ii) v6(ii, :) v5(ii)], 12, [t 'T{ii,:} (double)']);
    t_is(T2{:, :}, [v1 v2 v4 v6 v5], 12, [t 'T{:,:} (double)']);
    t_ok(isequal(T6{ii, 6:7}, ex(ii, :)), [t 'T{ii,j1:j2} (cell)']);
    t_ok(isequal(T6{:, 6:7}, ex), [t 'T{:,j1:j2} (cell)']);
    t_ok(isequal(T3{ii, :}, ex(ii, :)), [t 'T{ii,:} (cell)']);
    t_ok(isequal(T3{:, :}, ex), [t 'T{:,:} (cell)']);

    t = sprintf('%s : more subsasgn {} : ', cls);
    T2{ii, :} = [var1(ii) var2(ii) var4(ii) var6(ii, :) var5(ii)];
    t_is(T2{ii, :}, [var1(ii) var2(ii) var4(ii) var6(ii, :) var5(ii)], 12, [t 'T{ii,:} (double)']);
    T2{:, :} = [v1 v2 v4 v6 v5];
    t_is(T2{:, :}, [v1 v2 v4 v6 v5], 12, [t 'T{:,:} (cell)']);
    v1 = horzcat(var3, v3);
    v2 = horzcat(v3, var3);
    T6{ii, 6:7} = v1(ii, :);
    t_ok(isequal(T6{ii, 6:7}, v1(ii, :)), [t 'T{ii,j1:j2} (cell)']);
    T6{:, 6:7} = v2;
    t_ok(isequal(T6{:, 6:7}, v2), [t 'T{:,j1:j2} (cell)']);
    T3{ii, :} = v1(ii, :);
    t_ok(isequal(T3{ii, :}, v1(ii, :)), [t 'T{ii,:} (cell)']);
    T3{:, :} = v2;
    t_ok(isequal(T3{:, :}, v2), [t 'T{:,:} (cell)']);

    %% value class, not handle class
    t = sprintf('%s : value class, not handle', cls);
    T2 = T;
    T.igr(2) = 2;
    t_is(T.igr(2), 2, 12, t);
    t_is(T2.igr(2), 55, 12, t);

    %% nested tables
    for k2 = 1:nc
        nested_table_class = table_classes{k2};
        nstcls = upper(class_names{k2});

        t = sprintf('%s . %s : nested tables : ', cls, nstcls);
        T1 = nested_table_class([1;2;3],[4;5;6]);
        T2 = nested_table_class({'one';'two';'three'},{'four';'five';'six'});
        if have_feature('octave') && isa(T, 'mp_table_subclass') && ...
                                     isa(T1, 'mp_foo_table')
            t_skip(8, 'avoid https://savannah.gnu.org/bugs/index.php?65037');
        else
            T3 = table_class([10;20;30],{'uno';'dos';'tres'}, T1, T2, ...
                'VariableNames', {'v1', 'v2', 'T1', 'T2'});
            t_ok(isequal(T3.T1, T1), [t 'T.T1 == T1']);
            t_ok(isequal(T3.T2, T2), [t 'T.T2 == T2']);
            if have_feature('matlab') && isa(T1, 'mp_table_subclass') && ...
                                         isa(T3, 'mp_table_subclass')
                t_skip(4, 'multiple indexing of nested mp_table_subclasses');
            else
                t_ok(isequal(T3.T1([1;3], :), T1([1;3], :)), [t 'T.T1(ii, :) == T1(ii, :)']);
                t_ok(isequal(T3.T2([1;3], :), T2([1;3], :)), [t 'T.T2(ii, :) == T2(ii, :)']);
                t_ok(isequal(T3.T1.Var1, T1.Var1), [t 'T.T1.Var1 == T1.Var1']);
                t_ok(isequal(T3.T2.Var2([1;3], :), T2.Var2([1;3], :)), [t 'T.T2.Var2(ii) == T.T2.Var2(ii)']);
            end
            T4 = T3(:, :);
            T5 = T3([1;2;3], :);
            t_ok(isequal(T3, T4), [t 'T == T(:, :))']);
            t_ok(isequal(T3, T5), [t 'T == T(ii, :))']);
        end
    end
end

if nargout
    obj = T;
end

t_end;

function Tv = table_values(T)
if isa(T, 'mp_table_subclass')
    T = get_table(T);
end
if isa(T, 'table') && have_feature('matlab')
    Tv = cellfun(@(c)T{:, c}, num2cell(1:size(T, 2)), 'UniformOutput', false);
else
    Tv = T.Properties.VariableValues;
end

function show_me(T)
if isa(T, 'table') && have_feature('octave')
    prettyprint(T)
else
    T
end
