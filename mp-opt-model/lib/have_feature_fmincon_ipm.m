function [TorF, vstr, rdate] = have_feature_fmincon_ipm()
%HAVE_FEATURE_FMINCON_IPM  Detect availability/ver info for FMINCON w/Int Pt Mtd
%
%   Feature detection function implementing 'fmincon_ipm' tag for HAVE_FEATURE
%   to detect availability/version of FMINCON with Interior Point solver, from
%   the MATLAB Optimization Toolbox 4.x and later.
%
%   See also HAVE_FEATURE, NLPS_MASTER, FMINCON.

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

TorF = 0;
vstr = '';
rdate = '';
if have_feature('matlab')
    v = have_feature('fmincon', 'all');
    if v.av && v.vnum >= 4          %% Opt Tbx 4.0+ (R2008a+, MATLAB 7.6+)
        TorF = 1;
        vstr = v.vstr;
        rdate = v.date;
    end
end
