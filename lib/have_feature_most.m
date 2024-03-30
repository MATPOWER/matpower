function [TorF, vstr, rdate] = have_feature_most()
% have_feature_most - Detect availability/version info for MOST.
%
% Private feature detection function implementing ``'most'`` tag for
% have_feature to detect availability/version of MOST (MATPOWER Optimal Scheduling Tool).
%
% See also have_feature, most.

%   MATPOWER
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
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
