function [TorF, vstr, rdate] = have_feature_isequaln()
% have_feature_isequaln - Detect availability/version info for :func:`isequaln`.
%
% Private feature detection function implementing ``'isequaln'`` tag for
% have_feature to detect support for :func:`isequaln` function.
%
% See also have_feature, isequaln.

%   MP-Opt-Model
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
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
