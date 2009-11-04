function t_ok(cond, msg)
%T_OK  Tests if a condition is true.
%   t_ok(expr, msg) increments the global test count and if the expr
%   is true it increments the passed tests count, otherwise increments
%   the failed tests count. Prints 'ok' or 'not ok' followed by the
%   msg, unless the global variable t_quiet is true. Intended to be
%   called between calls to t_begin and t_end.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

global t_quiet;
global t_counter;
global t_ok_cnt;
global t_not_ok_cnt;

if nargin < 2 | strcmp(msg, '')
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
