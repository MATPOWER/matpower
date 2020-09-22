function [TorF, vstr, rdate] = have_feature_sdp_pf()
%HAVE_FEATURE_SDP_PF  Detect availability/version info for SDP_PF
%
%   Used by HAVE_FEATURE.

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

TorF = have_feature('yalmip') && exist('mpoption_info_sdp_pf', 'file') == 2;
if TorF
    v = sdp_pf_ver('all');
    vstr = v.Version;
    rdate = v.Date;
else
    vstr = '';
    rdate = '';
end
