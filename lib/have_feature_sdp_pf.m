function [TorF, vstr, rdate] = have_feature_sdp_pf()
% have_feature_sdp_pf - Detect availability/version info for SDP_PF.
%
% Private feature detection function implementing ``'sdp_pf'`` tag for
% have_feature to detect availability/version of SDP_PF, a |MATPOWER|
% extension for applications of semi-definite programming relaxations of
% power flow equations (https://github.com/MATPOWER/mx-sdp_pf/).
%
% See also have_feature.

%   MATPOWER
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
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
