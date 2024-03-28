function [TorF, vstr, rdate] = have_feature_knitro()
% have_feature_knitro - Detect availability/version info for Artelys Knitro.
%
% Private feature detection function implementing ``'knitro'`` tag for
% have_feature to detect availability/version of Artelys Knitro, a
% nonlinear programming solver (https://www.artelys.com/solvers/knitro/).
%
% See also have_feature, nlps_master, knitromatlab.

%   MP-Opt-Model
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

tmp = have_feature('knitromatlab', 'all');
if tmp.av
    TorF = tmp.av;
    vstr = tmp.vstr;
    rdate = tmp.date;
else
    tmp = have_feature('ktrlink', 'all');
    if tmp.av
        TorF = tmp.av;
        vstr = tmp.vstr;
        rdate = tmp.date;
    else
        TorF = 0;
        vstr = '';
        rdate = '';
    end
end
