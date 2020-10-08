function [TorF, vstr, rdate] = have_feature_quadprog_ls()
%HAVE_FEATURE_QUADPROG_LS  Detect availability/version info for QUADPROG w/large scale
%
%   Feature detection function implementing 'quadprog_ls' tag for HAVE_FEATURE
%   to detect availability/version of QUADPROG with support for the
%   large-scale interior point convex solver, from the MATLAB Optimization
%   Toolbox 6.x and later.
%
%   See also HAVE_FEATURE, HAVE_FEATURE_QUADPROG, QPS_MASTER, QUADPROG.

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
    v = have_feature('quadprog', 'all');
    if v.av && v.vnum >= 6          %% Opt Tbx 6.0+ (R2011a+, MATLAB 7.12+)
        TorF = 1;
        vstr = v.vstr;
        rdate = v.date;
    end
end
