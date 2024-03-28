function t_skip(cnt, msg)
% t_skip - Skips a number of tests.
% ::
%
%   t_skip(cnt, msg)
%
% Inputs:
%   cnt (integer) : number of tests to skip
%   msg (char array) : message to display for this set of skipped tests
%
% Increments the global test count and skipped tests count. Prints
% *"skipped x..y : "* followed by the ``msg``, unless t_begin was called
% with input ``quiet`` equal true.
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
% See also t_ok, t_is, t_file_match, t_str_match, t_begin, t_end, t_run_tests.

%   MP-Test
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
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
