function [TorF, vstr, rdate] = have_feature_most()
%HAVE_FEATURE_MOST  Detect availability/version info for MOST
%
%   Feature detection function implementing 'most' tag for HAVE_FEATURE
%   to detect availability/version of MOST (MATPOWER Optimal Scheduling Tool).
%
%   See also HAVE_FEATURE, MOST.

%   MATPOWER
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

TorF = exist('most', 'file') == 2;
if TorF
    v = mostver('all');
    vstr = v.Version;
    rdate = v.Date;
else
    vstr = '';
    rdate = '';
end
