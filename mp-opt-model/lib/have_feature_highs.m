function [TorF, vstr, rdate] = have_feature_highs()
% have_feature_highs - Detect availability/version info for HiGHS.
%
% Private feature detection function implementing ``'highs'`` tag for
% have_feature to detect availability/version of HiGHS optimizer
% (https://highs.dev).
%
% See also have_feature, qps_master, miqps_master, callhighs.

%   MP-Opt-Model
%   Copyright (c) 2004-2025, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

TorF = exist('callhighs', 'file') == 2 && exist('highsmex', 'file') == 3;
vstr = '';
rdate = '';
if TorF
    try
        vstr = char(callhighs(string('ver')));
    catch
        TorF = 0;
        fprintf('HiGHS Error!\n');
    end
end
