function t_test_fcns(quiet)
%T_TEST_FCNS  Test T_OK and T_IS and manually check output of failed tests.

%   MP-Test
%   Copyright (c) 2015, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

if nargin < 1
    quiet = 0;
end

ntests = 5;
npass = 9;
nfail = 12;

g = [];
e = [];

tol = 12;
k = 0;
str = {'EXPECT FAIL'; 'EXPECT PASS'};
f = @(s) str{s+1};

t_begin(npass+nfail, quiet);

k = k + 1; e(k) = 1;
t = sprintf('%s : t_ok(1, ...)', f(e(k)));
g(k) = t_ok(1, t);

k = k + 1; e(k) = 0;
t = sprintf('%s : t_ok(0, ...)', f(e(k)));
g(k) = t_ok(0, t);

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
