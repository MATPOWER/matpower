function [TorF, vstr, rdate] = have_feature_tralmopf()
% have_feature_tralmopf - Detect availability/version info for TRALMOPF
%
% Private feature detection function implementing ``'tralmopf'`` tag for
% have_feature to detect availability/version of TRALMOPF, trust region based
% augmented Langrangian optimal power flow (OPF) solver included in TSPOPF.
% (https://www.pserc.cornell.edu/tspopf)
%
% See also have_feature, tralmopf.

%   MATPOWER
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if have_feature('matlab') && exist('tralmopf', 'file') == 3;
    TorF = 1;
    v = tralmopfver('all');
    vstr = v.Version;
    rdate = v.Date;
else
    TorF = 0;
    vstr = '';
    rdate = '';
end
