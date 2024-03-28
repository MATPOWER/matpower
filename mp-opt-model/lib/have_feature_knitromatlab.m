function [TorF, vstr, rdate] = have_feature_knitromatlab()
% have_feature_knitromatlab - Detect availability/version info for :func:`knitromatlab`.
%
% Private feature detection function implementing ``'knitromatlab'`` tag for
% have_feature to detect availability/version of Artelys Knitro 9.0.0 and
% later with :func:`knitromatlab` or :func:`knitro_nlp` function.
%
% See also have_feature, have_feature_knitro, knitromatlab.

%   MP-Opt-Model
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% knitromatlab for Knitro 13.1 or greater
fcn_name = 'knitro_nlp';
TorF = exist(fcn_name, 'file');
if ~TorF
    %% knitromatlab for Knitro 9.0 - 12.0.x?
    fcn_name = 'knitromatlab';
    TorF = exist(fcn_name, 'file');
end
vstr = '';
rdate = '';
if TorF
    try
        str = evalc(sprintf('[x fval] = %s(@(x)1,1);', fcn_name));
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
