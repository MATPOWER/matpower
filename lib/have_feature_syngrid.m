function [TorF, vstr, rdate] = have_feature_syngrid()
% have_feature_syngrid - Detect availability/version info for SynGrid.
%
% Private feature detection function implementing ``'syngrid'`` tag for
% have_feature to detect availability/version of SynGrid, Synthetic Grid
% Creation for |MATPOWER| (https://github.com/MATPOWER/mx-syngrid).
%
% See also have_feature, syngrid.

%   MATPOWER
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

TorF = exist('syngrid', 'file') == 2;
if TorF
    v = sgver('all');
    vstr = v.Version;
    rdate = v.Date;
else
    vstr = '';
    rdate = '';
end
