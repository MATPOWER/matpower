function ok = t_is(got, expected, prec, msg)
% t_is - Tests if two matrices are identical to some tolerance.
% ::
%
%   t_is(got, expected, prec, msg)
%   ok = t_is(got, expected, prec, msg)
%
% Test passes if the maximum difference between corresponding elements of
% ``got`` and ``expected`` is less than :math:`10^{-prec}`.
%
% Inputs:
%   got (double) : numerical matrix of actual test results
%   expected (double) : numerical matrix of expected test results
%   prec (double) : defines the match tolerance, :math:`10^{-prec}`
%   msg (char array) : message to display for this test
%
% Output:
%   ok (boolean) : *(optional)* true if test passed, false if failed
%
% Increments the global test count and, if the test passes, then it
% increments the passed tests count, otherwise  increments the failed
% tests count. Prints *"ok"* or *"not ok"* followed by the ``msg``,
% unless t_begin was called with input ``quiet`` equal true.
%
% The input values can be real or complex, and they can be scalar, vector,
% or 2-d or higher matrices. If ``got`` is a vector or matrix and
% ``expected`` is a scalar or *NaN*, all elements must match the scalar.
% *NaN* values are considered to be equal to each other.
%
% Intended to be called between calls to t_begin and t_end.
%
% Example::
%
%   quiet = 0;
%   t_begin(5, quiet);
%   t_ok(pi > 3, 'size of pi');
%   t_skip(3, 'not yet written');
%   t_is(2+2, 4, 12, '2+2 still equals 4');
%   t_end;
%
% See also t_ok, t_file_match, t_str_match, t_skip, t_begin, t_end, t_run_tests.

%   MP-Test
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
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

if isreal(got) && isreal(expected)
    cplx = 0;
else
    cplx = 1;
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
            if cplx
                fprintf('\n%14s  %16s  %16s  %16g', ...
                    idxstr, format_complex(full(got(k(u))), '%g'), ...
                            format_complex(full(ex), '%g'), full(abs(got_minus_expected(k(u)))));
            else
                fprintf('\n%14s  %16g  %16g  %16g', ...
                    idxstr, full(got(k(u))), full(ex), full(abs(got_minus_expected(k(u)))));
            end
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

function s = format_complex(v, fmt)
if imag(v) < 0
    sgn = ' - ';
else
    sgn = ' + ';
end
s = sprintf([fmt sgn fmt 'i'], real(v), abs(imag(v)));
