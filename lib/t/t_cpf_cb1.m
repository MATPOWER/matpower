function [nx, cx, done, rollback, evnts, cb_data, results] = t_cpf_cb1(...
        k, nx, cx, px, done, rollback, evnts, cb_data, cb_args, results)
%T_CPF_CB1  User callback function 1 for continuation power flow testing.

%   MATPOWER
%   Copyright (c) 2013-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%-----  INITIAL call  -----
if k == 0
    cxx = struct(   'initial', 1, ...
                    'iteration', 0, ...
                    'final', 0  );
    nxx = cxx;
    cx.cb.cb1 = cxx;
    nx.cb.cb1 = nxx;
    if ~isfield(cx.cb, 'shared')
        cx.cb.shared = '';
        nx.cb.shared = '';
    end
    cx.cb.shared = [cx.cb.shared '1'];
    nx.cb.shared = [nx.cb.shared '1'];
else
    nxx = nx.cb.cb1;            %% get next callback state
    %%-----  ITERATION call  -----
    if k > 0
        nxx.iteration = nxx.iteration + 1;
        nx.cb.cb1 = nxx;        %% update next callback state
        nx.cb.shared = [nx.cb.shared '1'];  %% update next callback state
    %%-----  FINAL call  -----
    else    % k < 0
        results.cb1.initial     = nxx.initial;
        results.cb1.iteration   = nxx.iteration;
        results.cb1.final       = 1;
        results.shared          = nx.cb.shared;
    end
end