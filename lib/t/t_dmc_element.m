function t_dmc_element(quiet)
% t_dmc_element - Tests for mp.dmc_element.

%   MATPOWER
%   Copyright (c) 2022-2024, Power Systems Engineering Research Center (PSERC)
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

table_class = mp_table_class();

%%  1       2       3       4       5       6
%%  alpha   beta    uid     gamma   delta   epsilon
foo = [ ...
    11      5       1       0.1     9       0
    22      4       2       0.3     3       2
    33      3       3       0.5     7       0
    44      2       4       0.7     5       4
    55      1       5       0.9     1       6
];
baz = {
    'uno', 'one';
    'dos', 'two';
    'tres', 'three';
    'cuatro', 'four';
    'cinco', 'five';
};
d0 = struct('foo', foo, 'bar', struct('baz', {baz}));
ridx1 = 2;
ridx2 = [3;1];

r = [2;4;5];
var_names = {'uid', 'name', 'status', 'source_uid', ...
        'ids', 'twos', 'alpha', 'beta', 'gamma', 'delta', 'epsilon'};
ec = cell(5,1); [ec{:}] = deal('');
vals = {foo(:,3), baz(:,2), ones(5,1), ec, foo(:,3), 2*ones(5,1), ...
    foo(:,1), 10*foo(:,2), 2*foo(:,4), 5*foo(:,5), foo(:,6)};
tab0 = table_class(vals{:}, 'VariableNames', var_names);
tab1 = tab0(r,:); tab1.ids = [1:3]'; tab1.source_uid = r;

t_begin(24, quiet);

%% 1 - dme table rows correspond to original data rows
k = 'A';
t = sprintf('%s : dmce : ', k);
dmce = mp.dmce_test_widget();
t_ok(isa(dmce, 'mp.dmc_element'), [t 'class']);
t_is(dmce.use_r, 0, 12, [t 'dmce.use_r == 0']);

t = sprintf('%s : dme : ', k);
dme = mp.dme_test_widget();
t_ok(isa(dme, 'mp.dm_element'), [t 'class']);
t_ok(isempty(dme.tab), [t 'dme.tab empty']);

t = sprintf('%s : import(dme, d) : ', k);
dme = dmce.import(mp.dme_test_widget(), d0);
t_ok(~isempty(dme.tab), [t 'dme.tab not empty']);
t_ok(isequal(dme.tab, tab0), [t 'dme.tab']);

t = sprintf('%s : export(dme, d) : ', k);
d = struct();
d = dmce.export(dme, d);
ed = d0; ed.bar.baz(:,1) = ed.bar.baz(:,2); [ed.bar.baz{:,2}] = deal([]);
t_ok(isequal(d, ed), [t 'd']);

t = sprintf('%s : import(dme, d, var_names) : ', k);
d = d0; d.foo = d0.foo(:, end:-1:1);    %% reverse columns
tab = tab0; tmp = tab.alpha; tab.alpha = tab.epsilon; tab.epsilon = tmp;
dme = dmce.import(mp.dme_test_widget(), d0);
dme = dmce.import(dme, d, {'alpha', 'epsilon'});
t_ok(isequal(dme.tab, tab), [t 'dme.tab']);

t = sprintf('%s : export(dme, d, var_names) : ', k);
tmp = dme.tab.beta; dme.tab.beta = dme.tab.delta; dme.tab.delta = tmp;
ed = d0; ed.foo(:, [2 5]) = ed.foo(:, [5 2]);
ed.foo(:, 2) = ed.foo(:, 2)/2;
ed.foo(:, 5) = ed.foo(:, 5)*2;
d = dmce.export(dme, d0, {'beta', 'delta'});
t_ok(isequal(d, ed), [t 'd']);

t = sprintf('%s : import(dme, d, var_names, <scalar>) : ', k);
d = d0; d.foo = d0.foo(:, end:-1:1);    %% reverse columns
tab = tab0; tmp = tab.alpha(ridx1); tab.alpha(ridx1) = tab.epsilon(ridx1); tab.epsilon(ridx1) = tmp;
dme = dmce.import(mp.dme_test_widget(), d0);
dme = dmce.import(dme, d, {'alpha', 'epsilon'}, ridx1);
t_ok(isequal(dme.tab, tab), [t 'dme.tab']);

t = sprintf('%s : export(dme, d, var_names, <scalar>) : ', k);
tmp = dme.tab.beta; dme.tab.beta = dme.tab.delta; dme.tab.delta = tmp;
ed = d0; ed.foo(ridx1, [2 5]) = ed.foo(ridx1, [5 2]);
ed.foo(ridx1, 2) = ed.foo(ridx1, 2)/2;
ed.foo(ridx1, 5) = ed.foo(ridx1, 5)*2;
d = dmce.export(dme, d0, {'beta', 'delta'}, ridx1);
t_ok(isequal(d, ed), [t 'd']);

t = sprintf('%s : import(dme, d, var_names, <vector>) : ', k);
d = d0; d.foo = d0.foo(:, end:-1:1);    %% reverse columns
tab = tab0; tmp = tab.alpha(ridx2); tab.alpha(ridx2) = tab.epsilon(ridx2); tab.epsilon(ridx2) = tmp;
dme = dmce.import(mp.dme_test_widget(), d0);
dme = dmce.import(dme, d, {'alpha', 'epsilon'}, ridx2);
t_ok(isequal(dme.tab, tab), [t 'dme.tab']);

t = sprintf('%s : export(dme, d, var_names, <vector>) : ', k);
tmp = dme.tab.beta; dme.tab.beta = dme.tab.delta; dme.tab.delta = tmp;
ed = d0; ed.foo(ridx2, [2 5]) = ed.foo(ridx2, [5 2]);
ed.foo(ridx2, 2) = ed.foo(ridx2, 2)/2;
ed.foo(ridx2, 5) = ed.foo(ridx2, 5)*2;
d = dmce.export(dme, d0, {'beta', 'delta'}, ridx2);
t_ok(isequal(d, ed), [t 'd']);

%% 2 - dme table has r for indexing original data rows
k = 'B';
t = sprintf('%s : dmce : ', k);
dmce = mp.dmce_test_widget();
dmce.use_r = 1;
t_ok(isa(dmce, 'mp.dmc_element'), [t 'class']);
t_is(dmce.use_r, 1, 12, [t 'dmce.use_r == 1']);

t = sprintf('%s : import(dme, d) : ', k);
dme = dmce.import(mp.dme_test_widget(), d0);
t_ok(~isempty(dme.tab), [t 'dme.tab not empty']);
t_ok(isequal(dme.tab, tab1), [t 'dme.tab']);

t = sprintf('%s : export(dme, d) : ', k);
d = d0; d.foo = 0 * d.foo;
d = dmce.export(dme, d);
ed = d0;
ed.foo([1;3], :) = 0 * ed.foo([1;3], :);
ed.bar.baz(r,1) = ed.bar.baz(r,2);
t_ok(isequal(d, ed), [t 'd']);

t = sprintf('%s : import(dme, d, var_names) : ', k);
d = d0; d.foo(:, 1:end-2) = d0.foo(:, end-2:-1:1);  %% reverse columns
d.bar.baz(:, [1 2]) = d.bar.baz(:, [2 1]);          %% reverse columns
tab = tab1;
tmp = tab.alpha;
tab.alpha = tab.gamma / 2;
tab.gamma = tmp * 2;
tab.name = baz(r, 1);
dme = dmce.import(mp.dme_test_widget(), d0);
dme = dmce.import(dme, d, {'alpha', 'name', 'gamma'});
t_ok(isequal(dme.tab, tab), [t 'dme.tab']);

t = sprintf('%s : export(dme, d, var_names) : ', k);
tmp = dme.tab.beta; dme.tab.beta = dme.tab.delta; dme.tab.delta = tmp;
d = d0; d.foo = 0 * d.foo;
[d.bar.baz{:}] = deal('');
ed = d;
ed.foo(r, [2;5]) = d0.foo(r, [5;2]) * [0.5 0; 0 2];
ed.bar.baz(r,1) = d0.bar.baz(r,1);
d = dmce.export(dme, d, {'beta', 'delta', 'name'});
t_ok(isequal(d, ed), [t 'd']);

t = sprintf('%s : import(dme, d, var_names, <scalar>) : ', k);
d = d0; d.foo(:, 1:end-2) = d0.foo(:, end-2:-1:1);  %% reverse columns
d.bar.baz(:, [1 2]) = d.bar.baz(:, [2 1]);          %% reverse columns
tab = tab1;
tmp = tab.alpha(ridx1);
tab.alpha(ridx1) = tab.gamma(ridx1) / 2;
tab.gamma(ridx1) = tmp * 2;
tab.name(ridx1) = baz(r(ridx1), 1);
dme = dmce.import(mp.dme_test_widget(), d0);
dme = dmce.import(dme, d, {'alpha', 'name', 'gamma'}, ridx1);
t_ok(isequal(dme.tab, tab), [t 'dme.tab']);

t = sprintf('%s : export(dme, d, var_names, <scalar>) : ', k);
tmp = dme.tab.beta; dme.tab.beta = dme.tab.delta; dme.tab.delta = tmp;
d = d0; d.foo = 0 * d.foo;
[d.bar.baz{:}] = deal('');
ed = d;
ed.foo(r(ridx1), [2;5]) = d0.foo(r(ridx1), [5;2]) * [0.5 0; 0 2];
ed.bar.baz(r(ridx1),1) = d0.bar.baz(r(ridx1),1);
d = dmce.export(dme, d, {'beta', 'delta', 'name'}, ridx1);
t_ok(isequal(d, ed), [t 'd']);

t = sprintf('%s : import(dme, d, var_names, <vector>) : ', k);
d = d0; d.foo(:, 1:end-2) = d0.foo(:, end-2:-1:1);  %% reverse columns
d.bar.baz(:, [1 2]) = d.bar.baz(:, [2 1]);          %% reverse columns
tab = tab1;
tmp = tab.alpha(ridx2);
tab.alpha(ridx2) = tab.gamma(ridx2) / 2;
tab.gamma(ridx2) = tmp * 2;
tab.name(ridx2) = baz(r(ridx2), 1);
dme = dmce.import(mp.dme_test_widget(), d0);
dme = dmce.import(dme, d, {'alpha', 'name', 'gamma'}, ridx2);
t_ok(isequal(dme.tab, tab), [t 'dme.tab']);

t = sprintf('%s : export(dme, d, var_names, <vector>) : ', k);
tmp = dme.tab.beta; dme.tab.beta = dme.tab.delta; dme.tab.delta = tmp;
d = d0; d.foo = 0 * d.foo;
[d.bar.baz{:}] = deal('');
ed = d;
ed.foo(r(ridx2), [2;5]) = d0.foo(r(ridx2), [5;2]) * [0.5 0; 0 2];
ed.bar.baz(r(ridx2),1) = d0.bar.baz(r(ridx2),1);
d = dmce.export(dme, d, {'beta', 'delta', 'name'}, ridx2);
t_ok(isequal(d, ed), [t 'd']);

t_end;
