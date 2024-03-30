function [TorF, vstr, rdate] = have_feature_e4st()
% have_feature_e4st - Detect availability/version info for E4ST.
%
% Private feature detection function implementing ``'e4st'`` tag for
% have_feature to detect availability/version of E4ST, the Engineering,
% Economic, and Environmental Electricity Simulation Tool (https://e4st.com).
%
% See also have_feature.

%   MATPOWER
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

TorF = exist('e4st_ver', 'file') == 2;
if TorF
    v = e4st_ver('all');
    vstr = v.Version;
    rdate = v.Date;
else
    vstr = '';
    rdate = '';
end
