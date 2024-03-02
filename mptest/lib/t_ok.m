function ok = t_ok(cond, msg)
% t_ok - Tests whether a condition is true.
% ::
%
%   t_ok(expr, msg)
%   ok = t_ok(expr, msg)
%
% Test passes if ``expr`` evaluates to true.
%
% Inputs:
%   expr (boolean) : boolean expression from actual test results
%   msg (char array) : message to display for this test
%
% Output:
%   ok (boolean) : *(optional)* true if test passed, false if failed
%
% Increments the global test count and if the ``expr`` is true it
% increments the passed tests count, otherwise increments the failed
% tests count. Prints *"ok"* or *"not ok"* followed by the ``msg``,
% unless t_begin was called with input ``quiet`` equal true.
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
% See also t_is, t_file_match, t_str_match, t_skip, t_begin, t_end, t_run_tests.

%   MP-Test
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

global t_quiet;
global t_counter;
global t_ok_cnt;
global t_not_ok_cnt;

if nargin < 2 || strcmp(msg, '')
    msg = '';
else
    msg = [' - ', msg];
end
if cond
    t_ok_cnt = t_ok_cnt + 1;
else
    t_not_ok_cnt = t_not_ok_cnt + 1;
    if ~t_quiet
        fprintf('not ');
    end
end
if ~t_quiet
    fprintf('ok %d%s\n', t_counter, msg);
end
t_counter = t_counter + 1;
if nargout
    ok = cond;
end
