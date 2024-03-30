function [TorF, vstr, rdate] = have_feature_regexp_split()
% have_feature_regexp_split - Detect availability/version info for REGEXP ``'split'``.
%
% Private feature detection function implementing ``'regexp_split'`` tag for
% have_feature to detect support for the ``'split'`` argument to REGEXP.
%
% See also have_feature, regexp.

%   MATPOWER
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

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
