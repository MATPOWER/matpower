function t_test_fcns(quiet)
% t_test_fcns - Test MP-Test's individual test functions.
% ::
%
%   t_test_fcns
%   t_test_fcns(quiet)
%
% Test t_ok, t_is, t_str_match, t_file_match and manually check output of
% failed tests.

%   MP-Test
%   Copyright (c) 2015-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

if nargin < 1
    quiet = 0;
end

ntests = 5;
npass = 29;
nfail = 18;

g = [];
e = [];

tol = 12;
k = 0;
str = {'EXPECT FAIL'; 'EXPECT PASS'};
f = @(s) str{s+1};

t_begin(npass+nfail, quiet);

%% t_ok
k = k + 1; e(k) = 1;
t = sprintf('%s : t_ok(1, ...)', f(e(k)));
g(k) = t_ok(1, t);

k = k + 1; e(k) = 0;
t = sprintf('%s : t_ok(0, ...)', f(e(k)));
g(k) = t_ok(0, t);

%% t_is
k = k + 1; e(k) = 1;
t = sprintf('%s : t_is([], [], ...)', f(e(k)));
expected = [];
got = [];
g(k) = t_is(got, expected, tol, t);

k = k + 1; e(k) = 0;
t = sprintf('%s : t_is(1, [], ...)', f(e(k)));
expected = [];
got = 1;
g(k) = t_is(got, expected, tol, t);

k = k + 1; e(k) = 1;
t = sprintf('%s : t_is(3, 3, ...)', f(e(k)));
expected = 3;
got = 3;
g(k) = t_is(got, expected, tol, t);

k = k + 1; e(k) = 0;
t = sprintf('%s : t_is(1, 3, ...)', f(e(k)));
expected = 3;
got = 1;
g(k) = t_is(got, expected, tol, t);

k = k + 1; e(k) = 1;
t = sprintf('%s : t_is([3 3; 3 3], 3, ...)', f(e(k)));
expected = 3;
got = [3 3; 3 3];
g(k) = t_is(got, expected, tol, t);

k = k + 1; e(k) = 0;
t = sprintf('%s : t_is([3 4; 3 3], 3, ...)', f(e(k)));
expected = 3;
got = [3 4; 3 3];
g(k) = t_is(got, expected, tol, t);

k = k + 1; e(k) = 1;
t = sprintf('%s : t_is([3 3; 4 3], [3 3; 4 3], ...)', f(e(k)));
expected = [3 3; 4 3];
got = [3 3; 4 3];
g(k) = t_is(got, expected, tol, t);

k = k + 1; e(k) = 0;
t = sprintf('%s : t_is([3 4; 3 3], [3 3; 4 3], ...)', f(e(k)));
expected = [3 3; 4 3];
got = [3 4; 3 3];
g(k) = t_is(got, expected, tol, t);

k = k + 1; e(k) = 1;
t = sprintf('%s : t_is(ones(2,4,3), ones(2,4,3), ...)', f(e(k)));
expected = ones(2,4,3);
got = ones(2,4,3);
g(k) = t_is(got, expected, tol, t);

k = k + 1; e(k) = 0;
t = sprintf('%s : t_is(ones(2,4), ones(2,4,3), ...)', f(e(k)));
expected = ones(2,4,3);
got = ones(2,4);
g(k) = t_is(got, expected, tol, t);

k = k + 1; e(k) = 0;
t = sprintf('%s : t_is(ones(2,4,3), <one-elem-different> > tol, ...)', f(e(k)));
expected = ones(2,4,3);
got = ones(2,4,3);
got(2,3,2) = 1.0005;
g(k) = t_is(got, expected, 4, t);

k = k + 1; e(k) = 1;
t = sprintf('%s : t_is(ones(2,4,3), <one-elem-different> < tol, ...)', f(e(k)));
expected = ones(2,4,3);
got = ones(2,4,3);
got(2,3,2) = 1.0005;
g(k) = t_is(got, expected, 3, t);

k = k + 1; e(k) = 1;
t = sprintf('%s : t_is(NaN, NaN, ...)', f(e(k)));
expected = NaN;
got = NaN;
g(k) = t_is(got, expected, tol, t);

k = k + 1; e(k) = 0;
t = sprintf('%s : t_is([], NaN, ...)', f(e(k)));
expected = NaN;
got = [];
g(k) = t_is(got, expected, tol, t);

k = k + 1; e(k) = 0;
t = sprintf('%s : t_is(NaN, [], ...)', f(e(k)));
expected = [];
got = NaN;
g(k) = t_is(got, expected, tol, t);

k = k + 1; e(k) = 0;
t = sprintf('%s : t_is(1, NaN, ...)', f(e(k)));
expected = NaN;
got = 1;
g(k) = t_is(got, expected, tol, t);

k = k + 1; e(k) = 0;
t = sprintf('%s : t_is(NaN, 3, ...)', f(e(k)));
expected = 3;
got = NaN;
g(k) = t_is(got, expected, tol, t);

k = k + 1; e(k) = 1;
t = sprintf('%s : t_is(NaN(3,2), NaN, ...)', f(e(k)));
expected = NaN;
got = NaN(3,2);
g(k) = t_is(got, expected, tol, t);

k = k + 1; e(k) = 0;
t = sprintf('%s : t_is(int32(0.6), 1.0, ...)', f(e(k)));
expected = 0.6;
got = int32(1.0);
g(k) = t_is(got, expected, tol, t);

%% t_str_match
k = k + 1; e(k) = 1;
t = sprintf('%s : t_str_match(got, expected, ...)', f(e(k)));
expected = 'hello t_str_match';
got = 'hello t_str_match';
g(k) = t_str_match(got, expected, t);

k = k + 1; e(k) = 0;
t = sprintf('%s : t_str_match(got, expected, ...)', f(e(k)));
expected = 'hello t_str_match';
got = 'herro t_str_match';
g(k) = t_str_match(got, expected, t);

k = k + 1; e(k) = 1;
t = sprintf('%s : t_str_match(... reps[re])', f(e(k)));
expected = 'hello t_str_match';
got = 'herro t_str_match';
reps = {{'r{2}', 'll'}};
g(k) = t_str_match(got, expected, t, reps);

k = k + 1; e(k) = 0;
t = sprintf('%s : t_str_match(... reps[re])', f(e(k)));
expected = 'hello t_str_match';
got = 'herro t_str_match';
reps = {{'r{2}', 'LL'}};
g(k) = t_str_match(got, expected, t, reps);

k = k + 1; e(k) = 1;
t = sprintf('%s : t_str_match(... reps[!re])', f(e(k)));
expected = 'hello t_str_match';
got = 'herro t_str_match';
reps = {{'rr', 'll', 0}};
g(k) = t_str_match(got, expected, t, reps);

k = k + 1; e(k) = 0;
t = sprintf('%s : t_str_match(... reps[!re])', f(e(k)));
expected = 'hello t_str_match';
got = 'herro t_str_match';
reps = {{'r{2}', 'll', 0}};
g(k) = t_str_match(got, expected, t, reps);

k = k + 1; e(k) = 1;
t = sprintf('%s : t_str_match(... reps[both])', f(e(k)));
expected = 'hel4o t_str_match';
got = 'he56o t_str_match';
reps = {{'\d', 'l', 1, 1}};
g(k) = t_str_match(got, expected, t, reps);

k = k + 1; e(k) = 0;
t = sprintf('%s : t_str_match(... reps[both])', f(e(k)));
expected = 'heL4o t_str_match';
got = 'he56o t_str_match';
reps = {{'\d', 'l', 1, 1}};
g(k) = t_str_match(got, expected, t, reps);

%% t_file_match
[twd, n, ext] = fileparts(which('t_test_fcns'));
exp_fname = fullfile(twd, 't_have_feature', 'rithmaticker.m');
num = fix(1e7*rand);
fname1 = sprintf('test-output-1-%d.txt', num);
fname2 = sprintf('test-output-2-%d.txt', num);
fname3 = sprintf('test-output-3-%d.txt', num);
ex = fileread(exp_fname);
if isempty(strfind(ex, char([13 10])))
    ex2 = strrep(ex, char(10), char([13 10]));
else
    ex2 = strrep(ex, char([13 10]), char(10));
end
ex3 = regexprep(ex2, '\d', 'N');

k = k + 1; e(k) = 1;
t = sprintf('%s : fileread', f(e(k)));
g(k) = t_str_match(ex(1:8), 'function', t);

%% write files
filewrite(fname1, ex);
filewrite(fname2, ex2);
filewrite(fname3, ex3);

k = k + 1; e(k) = 1;
t = sprintf('%s : filewrite 1', f(e(k)));
ok = exist(fname1, 'file') == 2 && strcmp(ex, fileread(fname1));
g(k) = t_ok(ok, t);

k = k + 1; e(k) = 1;
t = sprintf('%s : filewrite 2', f(e(k)));
ok = exist(fname2, 'file') == 2 && ~strcmp(ex, fileread(fname2));
g(k) = t_ok(ok, t);

k = k + 1; e(k) = 1;
t = sprintf('%s : filewrite 3', f(e(k)));
ok = exist(fname3, 'file') == 2 && ~strcmp(ex, fileread(fname3));
g(k) = t_ok(ok, t);

%% t_file_match
k = k + 1; e(k) = 1;
t = sprintf('%s : t_file_match 1', f(e(k)));
g(k) = t_file_match(fname1, exp_fname, t);
k = k + 1; e(k) = 1;
g(k) = t_ok(exist(fname1, 'file') == 2, [t ' : file not deleted']);

k = k + 1; e(k) = 1;
t = sprintf('%s : t_file_match 1 (del)', f(e(k)));
g(k) = t_file_match(fname1, exp_fname, t, {}, 1);
k = k + 1; e(k) = 1;
g(k) = t_ok(exist(fname1, 'file') == 0, [t ' : file deleted']);

k = k + 1; e(k) = 1;
t = sprintf('%s : t_file_match 2', f(e(k)));
g(k) = t_file_match(fname2, exp_fname, t);
k = k + 1; e(k) = 1;
g(k) = t_ok(exist(fname2, 'file') == 2, [t ' : file not deleted']);

k = k + 1; e(k) = 1;
t = sprintf('%s : t_file_match 2 (del)', f(e(k)));
g(k) = t_file_match(fname2, exp_fname, t, {}, 1);
k = k + 1; e(k) = 1;
g(k) = t_ok(exist(fname2, 'file') == 0, [t ' : file deleted']);

k = k + 1; e(k) = 0;
t = sprintf('%s : t_file_match 3 (del)', f(e(k)));
g(k) = t_file_match(fname3, exp_fname, t, {}, 1);
k = k + 1; e(k) = 1;
g(k) = t_ok(exist(fname3, 'file') == 2, [t ' : file not deleted']);

k = k + 1; e(k) = 0;
t = sprintf('%s : t_file_match 3 (reps,del)', f(e(k)));
reps = {{'\d', 'N', 1}};
g(k) = t_file_match(fname3, exp_fname, t, reps);
k = k + 1; e(k) = 1;
g(k) = t_ok(exist(fname3, 'file') == 2, [t ' : file not deleted']);

k = k + 1; e(k) = 1;
t = sprintf('%s : t_file_match 3 (both,del)', f(e(k)));
reps = {{'[\dN]', 'M', 1, 1}};
g(k) = t_file_match(fname3, exp_fname, t, reps, 1);
k = k + 1; e(k) = 1;
g(k) = t_ok(exist(fname3, 'file') == 0, [t ' : file deleted']);



%%-----  reset tests, if this test passes, everything is as expected  -----
t_begin(ntests, quiet);

t = 'number of tests as expected';
t_is(length(g), npass+nfail, tol, t);
t_is(length(e), npass+nfail, tol, t);

t = 'verify test consistency';
t_is(sum( e), npass, tol, t);
t_is(sum(~e), nfail, tol, t);

t = 'passes and fails all as expected';
t_is(g, e, tol, t);

t_end;


function filewrite(fname, content)
[fid, msg] = fopen(fname, 'w');
if fid == -1
    error('%s\nt_test_fcns/filewrite: could not write file ''%s''', msg, fname);
end
fprintf(fid, '%s', content);
fclose(fid);
