function [nn, cc, cb_data, terminate, results] = t_cpf_cb1(...
        k, nn, cc, pp, rollback, critical, terminate, ...
        cb_data, cb_args, results)
%T_CPF_CB1  User callback function 1 for continuation power flow testing.

%   MATPOWER
%   Copyright (c) 2013-2016 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%-----  INITIAL call  -----
if k == 0
    cc.x.cb1.initial = 1;
    cc.x.cb1.iteration = 0;
    cc.x.cb1.final = 0;
%%-----  ITERATION call  -----
elseif k > 0
    nn.x.cb1.iteration = nn.x.cb1.iteration + 1;
%%-----  FINAL call  -----
else    % k < 0
    results.cb1.initial     = nn.x.cb1.initial;
    results.cb1.iteration   = nn.x.cb1.iteration;
    results.cb1.final       = 1;
end
