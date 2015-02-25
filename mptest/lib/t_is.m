function ok = t_is(got, expected, prec, msg)
%T_IS  Tests if two matrices are identical to some tolerance.
%   T_IS(GOT, EXPECTED, PREC, MSG) increments the global test count
%   and if the maximum difference between corresponding elements of
%   GOT and EXPECTED is less than 10^(-PREC) then it increments the
%   passed tests count, otherwise increments the failed tests count.
%   Prints 'ok' or 'not ok' followed by the MSG, unless the global
%   variable t_quiet is true. Intended to be called between calls to
%   T_BEGIN and T_END.
%
%   Optionally returns a true or false value indicating whether or
%   not the test succeeded. NaN's are considered to be equal to each
%   other.
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
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://matpower.org/ for more info.

global t_quiet;

if nargin < 4
    msg = '';
end
if nargin < 3 || isempty(prec)
    prec = 5;
end
[m, n] = size(expected);
if all(size(got) == [m, n]) || all([m, n] == [1 1])
    if m == 0 || n == 0
        condition = true;
    else
        %% check for NaNs!
        gNaN = find(isnan(got(:)));
        eNaN = find(isnan(expected(:)));
        if (~isscalar(expected) && ...
                (length(gNaN) ~= length(eNaN) || sum(gNaN-eNaN) ~= 0)) || ...
            (isscalar(expected) && ...
                    (( isnan(expected) && ~all(isnan(got))) || ...
                     (~isnan(expected) && any(isnan(got)))) )
            condition = false;
            max_diff = -1;
        elseif all(all(isnan(got))) && all(all(isnan(got)))
            condition = true;
        else
            got_minus_expected = got - expected;
            max_diff = max(max(abs(got_minus_expected)));
            condition = ( max_diff < 10^(-prec) );
        end
    end
else
    condition = false;
    max_diff = 0;
end

t_ok(condition, msg);
if ~condition && ~t_quiet
    if max_diff > 0
        [i, j, v] = find(~(abs(got_minus_expected) < 10^(-prec)));
        k = i+(j-1)*m;
        [vv, kk] = max(abs(got_minus_expected(k)));
        fprintf('  row     col          got             expected          got - exp\n');
        fprintf('-------  ------  ----------------  ----------------  ----------------');
        for u = 1:length(i)
            if isscalar(expected)
                ex = expected;
            else
                ex = expected(k(u));
            end
            fprintf('\n%6d  %6d  %16g  %16g  %16g', ...
                [i(u) j(u) got(k(u)) ex got_minus_expected(k(u))]');
            if u == kk
                fprintf('  *');
            end
        end
        fprintf('\nmax diff @ (%d,%d) = %g > allowed tol of %g\n\n', ...
            i(kk), j(kk), max_diff, 10^(-prec));
    elseif max_diff == -1
        fprintf('    mismatch in locations of NaNs\n');
    else
        fprintf('    dimension mismatch:\n');
        fprintf('             got: %d x %d\n', size(got));
        fprintf('        expected: %d x %d\n\n', size(expected));
    end
end
if nargout
    ok = condition;
end
