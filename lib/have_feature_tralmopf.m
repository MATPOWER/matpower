function [TorF, vstr, rdate] = have_feature_tralmopf()
%HAVE_FEATURE_TRALMOPF  Detect availability/version info for TRALMOPF
%
%   Used by HAVE_FEATURE.

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

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
