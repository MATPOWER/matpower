function [TorF, vstr, rdate] = have_feature_scpdipmopf()
% have_feature_scpdipmopf - Detect availability/version info for SCPDIPMOPF.
%
% Private feature detection function implementing ``'scpdipmopf'`` tag for
% have_feature to detect availability/version of SCPDIPMOPF, step-controlled
% primal-dual interior point method optimal power flow (OPF) solver
% included in TSPOPF. (https://www.pserc.cornell.edu/tspopf)
%
% See also have_feature, scpdipmopf.

%   MATPOWER
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
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
