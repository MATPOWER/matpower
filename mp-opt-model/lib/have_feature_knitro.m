function [TorF, vstr, rdate] = have_feature_knitro()
%HAVE_FEATURE_KNITRO  Detect availability/version info for Artelys Knitro
%
%   Feature detection function implementing 'knitro' tag for HAVE_FEATURE
%   to detect availability/version of Artelys Knitro, a nonlinear
%   programming solver (https://www.artelys.com/solvers/knitro/).
%
%   See also HAVE_FEATURE, NLPS_MASTER, KNITROMATLAB.

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
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
