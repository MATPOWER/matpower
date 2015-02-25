function [cb_state, results] = t_cpf_cb2(k, V_c, lam_c, V_p, lam_p, ...
										cb_data, cb_state, cb_args, results)
%T_CPF_CB2  User callback function 2 for continuation power flow testing.

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
