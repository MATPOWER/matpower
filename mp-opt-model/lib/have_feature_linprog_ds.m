function [TorF, vstr, rdate] = have_feature_linprog_ds()
% have_feature_linprog_ds - Detect availblty/ver info for :func:`linprog` w/dual simplex.
%
% Private feature detection function implementing ``'linprog_ds'`` tag for
% have_feature to detect availability/version of :func:`linprog` with
% support for the dual simplex method, from the MATLAB Optimization
% Toolbox 7.1 (R2014b) and later.
%
% See also have_feature, have_feature_linprog, qps_master, linprog.

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
    v = have_feature('linprog', 'all');
    if v.av && v.vnum >= 7.001      %% Opt Tbx 7.1+ (R2014b+, MATLAB 8.4+)
        TorF = 1;
        vstr = v.vstr;
        rdate = v.date;
    end
end
