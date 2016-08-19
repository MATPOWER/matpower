function [cb_state, nn, cc, cb_data, terminate, results] = t_cpf_cb1(...
        cont_steps, nn, cc, pp, rollback, critical, terminate, ...
        cb_data, cb_state, cb_args, results)
%T_CPF_CB1  User callback function 1 for continuation power flow testing.

%   MATPOWER
%   Copyright (c) 2013-2016 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

k = cont_steps;

%%-----  INITIAL call  -----
if k == 0
    cb_state.cb1.initial = 1;
    cb_state.cb1.iteration = 0;
    cb_state.cb1.final = 0;
%%-----  FINAL call  -----
elseif k < 0
    results.cb1.final = 1;
%%-----  ITERATION call  -----
else
    cb_state.cb1.iteration = cb_state.cb1.iteration + 1;
end
