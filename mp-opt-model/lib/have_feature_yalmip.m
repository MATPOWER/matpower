function [TorF, vstr, rdate] = have_feature_yalmip()
%HAVE_FEATURE_YALMIP  Detect availability/version info for YALMIP
%
%   Feature detection function implementing 'yalmip' tag for HAVE_FEATURE
%   to detect availability/version of YALMIP modeling platform
%   (https://yalmip.github.io).
%
%   See also HAVE_FEATURE, YALMIP.

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

TorF = exist('yalmip','file') == 2;
vstr = '';
rdate = '';
if TorF
    vstr = yalmip('version');
    if length(vstr) == 8
        yr = str2num(vstr(1:4));
        mo = str2num(vstr(5:6));
        dy = str2num(vstr(7:8));
        rdate = datestr([yr mo dy 0 0 0], 'dd-mmm-yyyy');
    end
end
