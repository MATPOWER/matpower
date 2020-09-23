function [TorF, vstr, rdate] = have_feature_smartmarket()
%HAVE_FEATURE_SMARTMARKET  Detect availability/version info for SMARTMARKET
%
%   Feature detection function implementing 'smartmarket' tag for HAVE_FEATURE
%   to detect availability/version of RUNMARKET and related files for running
%   an energy auction, found under smartmarket in MATPOWER Extras.
%   (https://github.com/MATPOWER/matpower-extras/).
%
%   See also HAVE_FEATURE, RUNMARKET.

%   MATPOWER
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

TorF = exist('runmarket', 'file') == 2;
if TorF
    v = mpver('all');
    vstr = v.Version;
    rdate = v.Date;
else
    vstr = '';
    rdate = '';
end
