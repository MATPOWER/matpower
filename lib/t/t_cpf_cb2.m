function [nx, cx, cb_data, done, results] = t_cpf_cb2(...
        k, nx, cx, px, rollback, critical, done, ...
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
if k == 0
    cxx = struct(   'initial', cb_args.cb2.initial, ...
                    'iteration', 0, ...
                    'final', 0  );
    nxx = cxx;
    cx.cb.cb2 = cxx;
    nx.cb.cb2 = nxx;
else
    nxx = nx.cb.cb2;            %% get next callback state
    %%-----  ITERATION call  -----
    if k > 0
        nxx.iteration = nxx.iteration + cb_args.cb2.iteration;
        nx.cb.cb2 = nxx;        %% update next callback state
    %%-----  FINAL call  -----
    else    % k < 0
        results.cb2.initial     = nxx.initial;
        results.cb2.iteration   = nxx.iteration;
        results.cb2.final       = cb_args.cb2.final;
    end
end