function t_skip(cnt, msg)
%T_SKIP  Skips a number of tests.
%   T_SKIP(CNT, MSG) increments the global test count and skipped tests
%   count. Prints 'skipped x..y : ' followed by the MSG, unless
%   T_BEGIN was called with input QUIET equal true. Intended to be
%   called between calls to T_BEGIN and T_END.
%
%   Example:
%       quiet = 0;
%       t_begin(5, quiet);
%       t_ok(pi > 3, 'size of pi');
%       t_skip(3, 'not yet written');
%       t_is(2+2, 4, 12, '2+2 still equals 4');
%       t_end;
%
%   See also T_OK, T_IS, T_BEGIN, T_END, T_RUN_TESTS.


%   MP-Test
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

global t_quiet;
global t_counter;
global t_skip_cnt;

if nargin < 2 || strcmp(msg, '')
    msg = '';
else
    msg = [' - ', msg];
end

t_skip_cnt = t_skip_cnt + cnt;
if ~t_quiet
    if cnt == 1
        fprintf('skipped %d%s\n', t_counter, msg);
    else
        fprintf('skipped %d..%d%s\n', t_counter, t_counter+cnt-1, msg);
    end
end
t_counter = t_counter + cnt;
