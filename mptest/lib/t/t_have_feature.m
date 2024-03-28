function t_have_feature(quiet)
% t_have_feature - Test have_feature function.
% ::
%
%   t_have_feature
%   t_have_feature(quiet)

%   MP-Test
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

if nargin < 1
    quiet = 0;
end

n_tests = 24;

t_begin(n_tests, quiet);

%% save cache
saved = have_feature('all', 'get_cache');

%% test set_cache/get_cache
t = 'set_cache/get_cache';
e = struct('foo', 2, 'bar', 3, 'baz', 'boz');
have_feature(e, 'set_cache');
t_ok(isequal(have_feature('all', 'get_cache'), e), t);
have_feature(64, 'set_cache');
t_ok(isequal(have_feature('all', 'get_cache'), 64), t);

%% clear cache
have_feature(struct(), 'set_cache');

%% Matlab vs. Octave
if exist('OCTAVE_VERSION', 'builtin') == 5
    t_ok(have_feature('octave'), 'Octave');
    t_ok(~have_feature('matlab'), 'not Matlab');
else
    t_ok(~have_feature('octave'), 'not Octave');
    t_ok(have_feature('matlab'), 'Matlab');
end

t = '<MPOM>/lib/t/t_have_feature must not be in path';
t_ok(exist('rithmaticker') ~= 2, t);

%% find path to this test file
cwd = pwd;
[p, n, e] = fileparts(which('t_have_feature'));

%% initially not available
t = 'have_feature(''rithmaticker'')';
t_ok(have_feature('rithmaticker') == 0, [t ' : not available']);

%% switch dir so it is available
cd(fullfile(p, 't_have_feature'));
cwd2 = pwd;

t = '<MPOM>/lib/t/t_have_feature must be in path';
t_ok(exist('rithmaticker') == 2, t);

%% still not available (cached)
t = 'have_feature(''rithmaticker'')';
t_ok(have_feature('rithmaticker') == 0, [t ' : still not available (cached)']);

%% clear cache, check again
have_feature('rithmaticker', 'clear_cache');
t_ok(have_feature('rithmaticker') == 1, [t ' : available (cache cleared)']);

t = 'have_feature(''rithmaticker'', 0)';
t_ok(have_feature('rithmaticker', 0) == 0, [t ' : disabled']);

t = 'have_feature(''rithmaticker'', 1)';
t_ok(have_feature('rithmaticker', 1) == 1, [t ' : enabled']);

t = 'have_feature(''rithmaticker'', false)';
t_ok(have_feature('rithmaticker', false) == 0, [t ' : disabled']);

t = 'have_feature(''rithmaticker'', true)';
t_ok(have_feature('rithmaticker', true) == 1, [t ' : enabled']);

cd(cwd);
t = 'successful switch to original working directory';
t_str_match(cwd, pwd, t);

t = 'have_feature(''rithmaticker'')';
t_ok(have_feature('rithmaticker') == 1, [t ' : still available (cached)']);

t = 'have_feature(''rithmaticker'', ''vstr'')';
t_str_match(have_feature('rithmaticker', 'vstr'), '3.1.4', t);

t = 'have_feature(''rithmaticker'', ''vnum'')';
t_is(have_feature('rithmaticker', 'vnum'), 3.001004, 12, t);

t = 'have_feature(''rithmaticker'', ''date'')';
t_str_match(have_feature('rithmaticker', 'date'), '30-May-2019', t);

t = 'have_feature(''rithmaticker'', ''all'')';
rv = have_feature('rithmaticker', 'all');
t_ok(isstruct(rv), [t ' : isstruct']);
t_is(rv.av, 1, 12, [t ' : av']);
t_str_match(rv.vstr, '3.1.4', [t ' : vstr']);
t_is(rv.vnum, 3.001004, 12, [t ' : vnum']);
t_str_match(rv.date, '30-May-2019', [t ' : date']);

%% clear cache, check again
t = 'have_feature(''rithmaticker'')';
have_feature('all', 'clear_cache');
t_ok(have_feature('rithmaticker') == 0, [t ' : not available (cache cleared)']);

%% restore cache
have_feature(saved, 'set_cache');

t_end;
