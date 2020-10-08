function [TorF, vstr, rdate] = have_feature_tralmopf()
%HAVE_FEATURE_TRALMOPF  Detect availability/version info for TRALMOPF
%
%   Feature detection function implementing 'tralmopf' tag for HAVE_FEATURE
%   to detect availability/version of TRALMOPF, trust region based
%   augmented Langrangian optimal power flow (OPF) solver included in TSPOPF.
%   (https://www.pserc.cornell.edu/tspopf)
%
%   See also HAVE_FEATURE, TRALMOPF.

%   MATPOWER
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
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
