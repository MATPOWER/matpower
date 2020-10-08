function [TorF, vstr, rdate] = have_feature_ktrlink()
%HAVE_FEATURE_KTRLINK  Detect availability/version info for KTRLINK
%
%   Feature detection function implementing 'ktrlink' tag for HAVE_FEATURE
%   to detect availability/version of Artelys Knitro prior to version 9.0.0,
%   which required the MATLAB Optimization Toolbox.
%
%   See also HAVE_FEATURE, HAVE_FEATURE_KNITRO, KTRLINK.

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% ktrlink for pre-Knitro 9.0, requires Optim Toolbox
TorF = exist('ktrlink', 'file');
vstr = '';
rdate = '';
if TorF
    try
        str = evalc(['[x fval] = ktrlink(@(x)1,1);']);
    end
    TorF = exist('fval', 'var') && fval == 1;
    if TorF
        pat = 'KNITRO ([^\s]+)\n|Knitro ([^\s]+)\n';
        [s,e,tE,m,t] = regexp(str, pat);
        if ~isempty(t)
            vstr = t{1}{1};
        end
    end
end
