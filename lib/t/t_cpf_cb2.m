function [nx, cx, done, rollback, evnts, cb_data, results] = t_cpf_cb2(...
        k, nx, cx, px, done, rollback, evnts, cb_data, cb_args, results)
%T_CPF_CB2  User callback function 2 for continuation power flow testing.

%   MATPOWER
%   Copyright (c) 2013-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%-----  INITIAL call  -----
if k == 0
    cxx = struct(   'initial', cb_args.initial, ...
                    'iteration', 0, ...
                    'final', 0  );
    nxx = cxx;
    cx.cb.cb2 = cxx;
    nx.cb.cb2 = nxx;
    if ~isfield(cx.cb, 'shared')
        cx.cb.shared = '';
        nx.cb.shared = '';
    end
    cx.cb.shared = [cx.cb.shared '2'];
    nx.cb.shared = [nx.cb.shared '2'];
else
    nxx = nx.cb.cb2;            %% get next callback state
    %%-----  ITERATION call  -----
    if k > 0
        nxx.iteration = nxx.iteration + cb_args.iteration;
        nx.cb.cb2 = nxx;        %% update next callback state
        nx.cb.shared = [nx.cb.shared '2'];  %% update next callback state
    %%-----  FINAL call  -----
    else    % k < 0
        results.cb2.initial     = nxx.initial;
        results.cb2.iteration   = nxx.iteration;
        results.cb2.final       = cb_args.final;
        results.shared          = nx.cb.shared;
    end
end