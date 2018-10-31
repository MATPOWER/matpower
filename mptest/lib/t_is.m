function ok = t_is(got, expected, prec, msg)
%T_IS  Tests if two matrices are identical to some tolerance.
%   T_IS(GOT, EXPECTED, PREC, MSG) increments the global test count
%   and if the maximum difference between corresponding elements of
%   GOT and EXPECTED is less than 10^(-PREC) then it increments the
%   passed tests count, otherwise increments the failed tests count.
%   Prints 'ok' or 'not ok' followed by the MSG, unless the global
%   variable t_quiet is true. The input values can be real or complex,
%   and they can be scalar, vector, or 2-d or higher matrices. If GOT
%   is a vector or matrix and EXPECTED is a scalar or NaN, all elements
%   must match the scalar. Intended to be called between calls to
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

%   MP-Test
%   Copyright (c) 2004-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

global t_quiet;

if nargin < 4
    msg = '';
end
if nargin < 3 || isempty(prec)
    prec = 5;
end

%% make sure we don't try to compare a double with an int
%% (difference can appear to be zero when it's actually not)
if isinteger(got) && ~isinteger(expected)
    got = double(got);
end
if ~isinteger(got) && isinteger(expected)
    expected = double(expected);
end

edims = size(expected);
gdims = size(got);
if (length(edims) == length(gdims) && all(edims == gdims)) || ...
        all(edims == 1) && ~any(gdims == 0)
    if all(edims == 0)
        condition = true;
    else
        %% check for NaNs!
        gNaN = find(isnan(got(:)));
        eNaN = find(isnan(expected(:)));
        if (~isscalar(expected) && ...
                (length(gNaN) ~= length(eNaN) || sum(gNaN-eNaN) ~= 0)) || ...
            (isscalar(expected) && ...
                    (( isnan(expected) && ~all(isnan(got(:)))) || ...
                     (~isnan(expected) && any(isnan(got(:))))) )
            condition = false;
            max_diff = -1;
        elseif all(isnan(expected(:))) && all(isnan(got(:)))
            condition = true;
        else
            got_minus_expected = got - expected;
            max_diff = max(abs(got_minus_expected(:)));
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
        k = find(~(abs(got_minus_expected(:)) < 10^(-prec)));
        [vv, kk] = max(abs(got_minus_expected(k)));
        fprintf('    index              got             expected      abs(got - exp)\n');
        fprintf('---------------  ----------------  ----------------  ----------------');
        for u = 1:length(k)
            if isscalar(expected)
                ex = expected;
            else
                ex = expected(k(u));
            end
            if isscalar(got)
                idxstr = '(1)';
            else
                idx = cell(1, length(gdims));
                [idx{:}] = ind2sub(gdims, k(u));
                idxstr = sprintf('%d,', idx{1:end-1});
                idxstr = sprintf('(%s%d)', idxstr, idx{end});
            end
            fprintf('\n%14s  %16g  %16g  %16g', ...
                idxstr, full(got(k(u))), full(ex), full(abs(got_minus_expected(k(u)))));
            if u == kk
                fprintf('  *');
                idxstrkk = idxstr;
            end
        end
        fprintf('\nmax diff @ %s = %g > allowed tol of %g\n\n', ...
            idxstrkk, full(max_diff), 10^(-prec));
    elseif max_diff == -1
        fprintf('    mismatch in locations of NaNs\n');
    else
        gidxstr = sprintf('%d x ', gdims(1:end-1));
        gidxstr = sprintf('%s%d', gidxstr, gdims(end));
        eidxstr = sprintf('%d x ', edims(1:end-1));
        eidxstr = sprintf('%s%d', eidxstr, edims(end));
        fprintf('    dimension mismatch:\n');
        fprintf('             got: %s\n', gidxstr);
        fprintf('        expected: %s\n\n', eidxstr);
    end
end
if nargout
    ok = condition;
end
