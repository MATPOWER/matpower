function t_have_fcn(quiet)
%T_HAVE_FCN  Tests for HAVE_FCN.

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

n_tests = 20;

t_begin(n_tests, quiet);

%% save cache
saved = have_fcn('all', 'get_cache');

%% test set_cache/get_cache
t = 'set_cache/get_cache';
e = struct('foo', 2, 'bar', 3, 'baz', 'boz');
have_fcn(e, 'set_cache');
t_ok(isequal(have_fcn('all', 'get_cache'), e), t);
have_fcn(64, 'set_cache');
t_ok(isequal(have_fcn('all', 'get_cache'), 64), t);

%% clear cache
have_fcn(struct(), 'set_cache');

%% Matlab vs. Octave
if exist('OCTAVE_VERSION', 'builtin') == 5
    t_ok(have_fcn('octave'), 'Octave');
    t_ok(~have_fcn('matlab'), 'not Matlab');
else
    t_ok(~have_fcn('octave'), 'not Octave');
    t_ok(have_fcn('matlab'), 'Matlab');
end

t = '$MATPOWER/t/t_feval_w_path must not be in path';
t_ok(exist('rithmaticker') ~= 2, t);

%% find path to this test file
cwd = pwd;
[p, n, e] = fileparts(which('t_have_fcn'));

%% initially not available
t = 'have_fcn(''rithmaticker'')';
t_ok(have_fcn('rithmaticker') == 0, [t ' : not available']);

%% switch dir so it is available
cd(fullfile(p, 't_feval_w_path'));
cwd2 = pwd;

t = '$MATPOWER/t/t_feval_w_path must be in path';
t_ok(exist('rithmaticker') == 2, t);

%% still not available (cached)
t = 'have_fcn(''rithmaticker'')';
t_ok(have_fcn('rithmaticker') == 0, [t ' : still not available (cached)']);

%% clear cache, check again
have_fcn('rithmaticker', 'clear_cache');
t_ok(have_fcn('rithmaticker') == 1, [t ' : available (cache cleared)']);

cd(cwd);
t = 'successful switch to original working directory';
t_ok(strcmp(cwd, pwd), t);

t = 'have_fcn(''rithmaticker'')';
t_ok(have_fcn('rithmaticker') == 1, [t ' : still available (cached)']);

t = 'have_fcn(''rithmaticker'', ''vstr'')';
t_ok(strcmp(have_fcn('rithmaticker', 'vstr'), '3.1.4'), t);

t = 'have_fcn(''rithmaticker'', ''vnum'')';
t_is(have_fcn('rithmaticker', 'vnum'), 3.001004, 12, t);

t = 'have_fcn(''rithmaticker'', ''date'')';
t_ok(strcmp(have_fcn('rithmaticker', 'date'), '30-May-2019'), t);

t = 'have_fcn(''rithmaticker'', ''all'')';
rv = have_fcn('rithmaticker', 'all');
t_ok(isstruct(rv), [t ' : isstruct']);
t_is(rv.av, 1, 12, [t ' : av']);
t_ok(strcmp(rv.vstr, '3.1.4'), [t ' : vstr']);
t_is(rv.vnum, 3.001004, 12, [t ' : vnum']);
t_ok(strcmp(rv.date, '30-May-2019'), [t ' : date']);

%% clear cache, check again
t = 'have_fcn(''rithmaticker'')';
have_fcn('all', 'clear_cache');
t_ok(have_fcn('rithmaticker') == 0, [t ' : not available (cache cleared)']);

%% restore cache
have_fcn(saved, 'set_cache');

t_end;
