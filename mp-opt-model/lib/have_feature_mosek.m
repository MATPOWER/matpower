function [TorF, vstr, rdate] = have_feature_mosek()
% have_feature_mosek - Detect availability/version info for MOSEK.
%
% Private feature detection function implementing ``'mosek'`` tag for
% have_feature to detect availability/version of MOSEK, LP/QP/MILP/MIQP solver
% (https://www.mosek.com/).
%
% See also have_feature, qps_master, miqps_master, mosekopt.

%   MP-Opt-Model
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

TorF = exist('mosekopt', 'file') == 3;
vstr = '';
rdate = '';
if TorF
    % MOSEK Version 6.0.0.93 (Build date: 2010-10-26 13:03:27)
    % MOSEK Version 6.0.0.106 (Build date: 2011-3-17 10:46:54)
    % MOSEK Version 7.0.0.134 (Build date: 2014-10-2 11:10:02)
    pat = 'Version (\.*\d)+.*Build date: (\d+-\d+-\d+)';
    [s,e,tE,m,t] = regexp(evalc('mosekopt'), pat);
    if isempty(t)
        [r, res] = mosekopt('version');
        v = res.version;
        vstr = sprintf('%d.%d.%d.%d', ...
            v.major, v.minor, v.build, v.revision);
    else
        vstr = t{1}{1};
        rdate = datestr(t{1}{2}, 'dd-mmm-yyyy');
    end
end
