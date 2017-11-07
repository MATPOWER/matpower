function t_feval_w_path(quiet)
%T_FEVAL_W_PATH  Tests for FEVAL_W_PATH.

%   MATPOWER
%   Copyright (c) 2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

n_tests = 23;

t_begin(n_tests, quiet);

t = '$MATPOWER/t/t_feval_w_path must not be in path';
t_ok(exist('rithmaticker') ~= 2, t);

%% find path to this test file
cwd = pwd;
[p, n, e] = fileparts(which('t_feval_w_path'));

%% test with empty path (do NOT cd anywhere)
t = 'ab = feval_w_path(<empty path>, fname, a, b)';
rv = feval_w_path('', 't_ok', 1, 'true');
t_ok(strcmp(cwd, pwd), [t ' : cwd unchanged']);
t_is(rv, 1, 15, t);

%% test with absolute path
t = 'ab = feval_w_path(<absolute path>, fname, a, b)';
rv = feval_w_path(fullfile(p, 't_feval_w_path'), 'rithmaticker', 3, 4);
t_ok(strcmp(cwd, pwd), [t ' : cwd unchanged']);
t_is(rv, 12, 15, t);

t = 'abc = feval_w_path(<absolute path>, fname, a, b, c)';
rv = feval_w_path(fullfile(p, 't_feval_w_path'), 'rithmaticker', 3, 4, 5);
t_ok(strcmp(cwd, pwd), [t ' : cwd unchanged']);
t_is(rv, 60, 15, t);

t = 'abcd = feval_w_path(<absolute path>, fname, a, b, c, d)';
rv = feval_w_path(fullfile(p, 't_feval_w_path'), 'rithmaticker', 3, 4, 5, 6);
t_ok(strcmp(cwd, pwd), [t ' : cwd unchanged']);
t_is(rv, 360, 15, t);

t = '[ab, a_b] = feval_w_path(<absolute path>, fname, a, b)';
[rv, rv2] = feval_w_path(fullfile(p, 't_feval_w_path'), 'rithmaticker', 3, 4);
t_ok(strcmp(cwd, pwd), [t ' : cwd unchanged']);
t_is(rv, 12, 15, t);
t_is(rv2, 7, 15, t);

%% switch dir to be able to use relative path
cd(p);
cwd2 = pwd;

%% test with relative path
t = 'ab = feval_w_path(<relative path>, fname, a, b)';
rv = feval_w_path('t_feval_w_path', 'rithmaticker', 3, 4);
t_ok(strcmp(cwd2, pwd), [t ' : cwd unchanged']);
t_is(rv, 12, 15, t);

t = 'abc = feval_w_path(<relative path>, fname, a, b, c)';
rv = feval_w_path('t_feval_w_path', 'rithmaticker', 3, 4, 5);
t_ok(strcmp(cwd2, pwd), [t ' : cwd unchanged']);
t_is(rv, 60, 15, t);

t = 'abcd = feval_w_path(<relative path>, fname, a, b, c, d)';
rv = feval_w_path('t_feval_w_path', 'rithmaticker', 3, 4, 5, 6);
t_ok(strcmp(cwd2, pwd), [t ' : cwd unchanged']);
t_is(rv, 360, 15, t);

t = '[ab, a_b] = feval_w_path(<relative path>, fname, a, b)';
[rv, rv2] = feval_w_path('t_feval_w_path', 'rithmaticker', 3, 4);
t_ok(strcmp(cwd2, pwd), [t ' : cwd unchanged']);
t_is(rv, 12, 15, t);
t_is(rv2, 7, 15, t);

cd(cwd);
t = 'successful switch to original working directory';
t_ok(strcmp(cwd, pwd), t);

t_end;
