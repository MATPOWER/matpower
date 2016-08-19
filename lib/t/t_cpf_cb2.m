function [nn, cc, cb_data, terminate, results] = t_cpf_cb2(...
        cont_steps, nn, cc, pp, rollback, critical, terminate, ...
        cb_data, cb_args, results)
%T_CPF_CB2  User callback function 2 for continuation power flow testing.

%   MATPOWER
%   Copyright (c) 2013-2016 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%-----  INITIAL call  -----
if cont_steps == 0
    cc.x.cb2.initial = cb_args.cb2.initial;
    cc.x.cb2.iteration = 0;
    cc.x.cb2.final = 0;
%%-----  ITERATION call  -----
elseif cont_steps > 0
    nn.x.cb2.iteration = nn.x.cb2.iteration + cb_args.cb2.iteration;
%%-----  FINAL call  -----
else    % cont_steps < 0
    results.cb2.initial     = nn.x.cb2.initial;
    results.cb2.iteration   = nn.x.cb2.iteration;
    results.cb2.final       = cb_args.cb2.final;
end
