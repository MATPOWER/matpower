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
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

global t_quiet;

if nargin < 4
    msg = '';
end
if nargin < 3 || isempty(prec)
    prec = 5;
end
[m, n] = size(expected);
if all(size(got) == [m, n]) || all([m, n] == [1 1])
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
        [i, j, v] = find(abs(got_minus_expected) >= 10^(-prec));
        k = i+(j-1)*m;
        [vv, kk] = max(abs(got_minus_expected(k)));
        fprintf('  row     col          got             expected          got - exp\n');
        fprintf('-------  ------  ----------------  ----------------  ----------------');
        for u = 1:length(i)
            fprintf('\n%6d  %6d  %16g  %16g  %16g', ...
                [i(u) j(u) got(k(u)) expected(k(u)) got_minus_expected(k(u))]');
            if u == kk
                fprintf('  *');
            end
        end
        fprintf('\nmax diff @ (%d,%d) = %g > allowed tol of %g\n\n', ...
            i(kk), j(kk), max_diff, 10^(-prec));
    else
        fprintf('    dimension mismatch:\n');
        fprintf('             got: %d x %d\n', size(got));
        fprintf('        expected: %d x %d\n\n', size(expected));
    end
end
