function [TorF, vstr, rdate] = have_feature_fmincon_ipm()
% have_feature_fmincon_ipm - Detect availability/ver info for :func:`fmincon` w/Int Pt Mtd.
%
% Private feature detection function implementing ``'fmincon_ipm'`` tag for
% have_feature to detect availability/version of :func:`fmincon` with
% Interior Point solver, from the MATLAB Optimization Toolbox 4.x and later.
%
% See also have_feature, nlps_master, fmincon.

%   MP-Opt-Model
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
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
