function [cb_state, results] = t_cpf_cb2(k, V_c, lam_c, V_p, lam_p, ...
										cb_data, cb_state, cb_args, results)
%T_CPF_CB2  User callback function 2 for continuation power flow testing.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2013 by Power System Engineering Research Center (PSERC)
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

%%-----  INITIAL call  -----
if k == 0
	cb_state.cb2.initial = cb_args.cb2.initial;
	cb_state.cb2.iteration = 0;
	cb_state.cb2.final = 0;
%%-----  FINAL call  -----
elseif nargout == 2
	results.cb2.final = cb_args.cb2.final;
%%-----  ITERATION call  -----
else
	cb_state.cb2.iteration = cb_state.cb2.iteration + cb_args.cb2.iteration;
end
