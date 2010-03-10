function t_is(got, expected, prec, msg)
%T_IS  Tests if two matrices are identical to some tolerance.
%   T_IS(GOT, EXPECTED, PREC, MSG) increments the global test count
%   and if the maximum difference between corresponding elements of
%   GOT and EXPECTED is less than 10^(-PREC) then it increments the
%   passed tests count, otherwise increments the failed tests count.
%   Prints 'ok' or 'not ok' followed by the MSG, unless the global
%   variable t_quiet is true. Intended to be called between calls to
%   T_BEGIN and T_END.
%
%   Example:
%       quiet = 0;
%       t_begin(5, quiet);
%       t_ok(pi > 3, 'size of pi');
%       t_skip(3, 'not yet written');
%       t_is(2+2, 4, 12, '2+2 still equals 4');
%       t_end;
%
%   See also T_OK, T_SKIP, T_BEGIN, T_END, T_RUN_TESTS.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2004-2010 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

global t_quiet;

if nargin < 4
    msg = '';
end
if nargin < 3 || isempty(prec)
    prec = 5;
end

if all(size(got) == size(expected)) || all(size(expected) == [1 1])
    got_minus_expected = got - expected;
    max_diff = max(max(abs(got_minus_expected)));
    condition = ( max_diff < 10^(-prec) );
else
    condition = false;
    max_diff = 0;
end

t_ok(condition, msg);
if ~condition && ~t_quiet
    if max_diff ~= 0
        got
        expected
        got_minus_expected
        fprintf('max diff = %g (allowed tol = %g)\n\n', max_diff, 10^(-prec));
    else
        fprintf('    dimension mismatch:\n');
        fprintf('             got: %d x %d\n', size(got));
        fprintf('        expected: %d x %d\n\n', size(expected));
    end
end
