function t_skip(cnt, msg)
%T_SKIP  Skips a number of tests.
%   T_SKIP(CNT, MSG) increments the global test count and skipped tests
%   count. Prints 'skipped tests x..y : ' followed by the MSG, unless the
%   global variable t_quiet is true. Intended to be called between calls to
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
%   See also T_OK, T_IS, T_BEGIN, T_END, T_RUN_TESTS.


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
