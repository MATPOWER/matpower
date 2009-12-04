function t_skip(cnt, msg)
%T_SKIP  Skips a number of tests.
%   t_skip(cnt, msg) increments the global test count and skipped tests
%   count. Prints 'skipped tests x..y : ' followed by the msg, unless the
%   global variable t_quiet is true. Intended to be called between calls to
%   t_begin and t_end.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

global t_quiet;
global t_counter;
global t_skip_cnt;

if nargin < 2 || strcmp(msg, '')
    msg = '';
else
    msg = [' : ', msg];
end

t_skip_cnt = t_skip_cnt + cnt;
if ~t_quiet
    fprintf('skipped tests %d..%d%s\n', t_counter, t_counter+cnt-1, msg);
end
t_counter = t_counter + cnt;
