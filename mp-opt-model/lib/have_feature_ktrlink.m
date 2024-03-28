function [TorF, vstr, rdate] = have_feature_ktrlink()
% have_feature_ktrlink - Detect availability/version info for :func:`ktrlink`.
%
% Private feature detection function implementing ``'ktrlink'`` tag for
% have_feature to detect availability/version of Artelys Knitro prior to
% version 9.0.0, with the :func:`ktrlink` function, which required the
% MATLAB Optimization Toolbox.
%
% See also have_feature, have_feature_knitro, ktrlink.

%   MP-Opt-Model
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
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
