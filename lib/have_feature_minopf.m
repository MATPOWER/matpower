function [TorF, vstr, rdate] = have_feature_minopf()
% have_feature_minopf - Detect availability/version info for MINOPF.
%
% Private feature detection function implementing ``'minopf'`` tag for
% have_feature to detect availability/version of MINOPF, a MINOS-based
% optimal power flow (OPF) solver.
%
% See also have_feature, minopf.

%   MATPOWER
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

TorF = exist('minopf', 'file') == 3;
if TorF
    v = minopfver('all');
    vstr = v.Version;
    rdate = v.Date;
else
    vstr = '';
    rdate = '';
end
