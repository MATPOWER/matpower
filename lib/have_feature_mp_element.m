function [TorF, vstr, rdate] = have_feature_mp_element()
%HAVE_FEATURE_MP_ELEMENT  Detect whether to use MP-Element by default
%
%   Feature detection function implementing 'mp_element' tag for HAVE_FEATURE
%   to availability/version of MP-Element.
%
%   See also HAVE_FEATURE, MP_ELEMENT.

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

TorF = exist('mp_element') ~= 0;
v = mpver('all');
vstr = v.Version;
rdate = v.Date;
