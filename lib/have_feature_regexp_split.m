function [TorF, vstr, rdate] = have_feature_regexp_split()
%HAVE_FEATURE_REGEXP_SPLIT  Detect availability/version info for REGEXP 'split'
%
%   Used by HAVE_FEATURE.

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
    v = have_feature('matlab', 'all');
    if v.av
        TorF = 1;
        vstr = v.vstr;
        rdate = v.date;
    end
else
    v = have_feature('octave', 'all');
    if v.av
        TorF = 1;
        vstr = v.vstr;
        rdate = v.date;
    end
end
