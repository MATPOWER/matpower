function [TorF, vstr, rdate] = have_feature_highs_hipo()
% have_feature_highs_hipo - Detect availability/version info for HiPO solver for HiGHS.
%
% Private feature detection function implementing ``'highs_hipo'`` tag for
% have_feature to detect availability/version of HiPO solver within HiGHS
% optimizer (https://highs.dev).
%
% See also have_feature, qps_master, callhighs.

%   MP-Opt-Model
%   Copyright (c) 2004-2025, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.


v = have_feature('highs', 'all');
TorF = v.av;
vstr = v.vstr;
rdate = v.date;
if TorF
    try
        s = warning('query');
        warning('off');
        x = feval('qps_highs', [],[1; 1],[1 1],[2],[2],[1; 1],[1; 1],[1; 1],struct('verbose', 0,'highs_opt',struct('solver','hipo')));
    catch
        TorF = 0;
        warning(s);
    end
end
