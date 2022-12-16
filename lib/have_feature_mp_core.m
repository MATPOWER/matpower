function [TorF, vstr, rdate] = have_feature_mp_core()
%HAVE_FEATURE_MP_CORE  Detect availability of MP-Core
%
%   Feature detection function implementing 'mp_core' tag for HAVE_FEATURE
%   to detect availability/version of MP-Core.
%
%   See also HAVE_FEATURE.

%   MATPOWER
%   Copyright (c) 2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% requires MATLAB 9.1+ or Octave 6.2+
TorF = exist('+mp/nm_element.m') ~= 0 && ...
    ((have_feature('matlab') && have_feature('matlab', 'vnum') >= 9.001) || ...
     (have_feature('octave') && have_feature('octave', 'vnum') >= 6.002));
v = mpver('all');
vstr = v.Version;
rdate = v.Date;
