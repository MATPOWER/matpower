function [cb_state, results] = t_cpf_cb1(k, V_c, lam_c, V_p, lam_p, ...
										cb_data, cb_state, cb_args, results)
%T_CPF_CB1  User callback function 1 for continuation power flow testing.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2013 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://matpower.org/ for more info.


%%-----  INITIAL call  -----
if k == 0
	cb_state.cb1.initial = 1;
	cb_state.cb1.iteration = 0;
	cb_state.cb1.final = 0;
%%-----  FINAL call  -----
elseif nargout == 2
	results.cb1.final = 1;
%%-----  ITERATION call  -----
else
	cb_state.cb1.iteration = cb_state.cb1.iteration + 1;
end
