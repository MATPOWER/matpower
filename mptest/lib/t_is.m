function t_is(got, expected, prec, msg)
%T_IS  Tests if two matrices are identical to some tolerance.
%   t_is(got, expected, prec, msg) increments the global test count
%   and if the maximum difference between corresponding elements of
%   got and expected is less than 10^(-prec) then it increments the
%   passed tests count, otherwise increments the failed tests count.
%   Prints 'ok' or 'not ok' followed by the msg, unless the global
%   variable t_quiet is true. Intended to be called between calls to
%   t_begin and t_end.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 4
    msg = '';
end
if nargin < 3 | isempty(prec)
    prec = 5;
end

got_minus_expected = got - expected;
max_diff = max(max(abs(got_minus_expected)));
cond = ( max_diff < 10^(-prec) );

t_ok(cond, msg);
if ~cond
    got
    expected
    got_minus_expected
    fprintf('max diff = %g (allowed tol = %g)\n\n', max_diff, 10^(-prec));
end

return;
