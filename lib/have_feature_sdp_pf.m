function [TorF, vstr, rdate] = have_feature_sdp_pf()
%HAVE_FEATURE_SDP_PF  Detect availability/version info for SDP_PF
%
%   Feature detection function implementing 'sdp_pf' tag for HAVE_FEATURE
%   to detect availability/version of SDP_PF, a MATPOWER extension for
%   applications of semi-definite programming relaxations of power flow
%   equations (https://github.com/MATPOWER/mx-sdp_pf/).
%
%   See also HAVE_FEATURE.

%   MATPOWER
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

TorF = have_feature('yalmip') && exist('mpoption_info_sdp_pf', 'file') == 2;
if TorF
    v = sdp_pf_ver('all');
    vstr = v.Version;
    rdate = v.Date;
else
    vstr = '';
    rdate = '';
end
