function [TorF, vstr, rdate] = have_feature_isequaln()
%HAVE_FEATURE_ISEQUALN  Detect availability/version info for ISEQUALN
%
%   Feature detection function implementing 'isequaln' tag for HAVE_FEATURE
%   to detect support for ISEQUALN function.
%
%   See also HAVE_FEATURE, ISEQUALN.

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

TorF = exist('isequaln') ~= 0;
if TorF
    if have_feature('matlab')
        %% introduced after 7.10 (R2010a) and by 7.14 (R2012a)
        v = have_feature('matlab', 'all');
    else
        %% introduced in Octave 4.4 (not in 4.2.2)
        v = have_feature('octave', 'all');
    end
    vstr  = v.vstr;
    rdate = v.date;
else
    vstr = '';
    rdate = '';
end
