function [TorF, vstr, rdate] = have_feature_bpmpd()
%HAVE_FEATURE_BPMPD  Detect availability/version info for BPMPD
%
%   Feature detection function implementing 'bpmpd' tag for HAVE_FEATURE
%   to detect availability/version of BPMPD_MEX (interior point LP/QP solver).
%
%   See also HAVE_FEATURE, QPS_MASTER, BP, BPOPT.

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

TorF = exist('bp', 'file') == 3;
if TorF
    v = bpver('all');
    vstr = v.Version;
    rdate = v.Date;
else
    vstr = '';
    rdate = '';
end
