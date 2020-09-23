function [TorF, vstr, rdate] = have_feature_scpdipmopf()
%HAVE_FEATURE_SCPDIPMOPF  Detect availability/version info for SCPDIPMOPF
%
%   Feature detection function implementing 'scpdipmopf' tag for HAVE_FEATURE
%   to detect availability/version of SCPDIPMOPF, a step-controlled
%   primal-dual interior point method Optimal Power Flow solver.
%
%   See also HAVE_FEATURE, SCPDIPMOPF.

%   MATPOWER
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if have_feature('matlab') && exist('scpdipmopf', 'file') == 3;
    TorF = 1;
    v = scpdipmopfver('all');
    vstr = v.Version;
    rdate = v.Date;
else
    TorF = 0;
    vstr = '';
    rdate = '';
end
