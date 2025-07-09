function [TorF, vstr, rdate] = have_feature_mips()
% have_feature_mips - Detect availability/version info for MIPS.
%
% Private feature detection function implementing ``mips`` tag for
% have_feature to detect support for ``mips()`` function.
%
%   See also have_feature.

%   MIPS
%   Copyright (c) 2004-2025, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MIPS.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mips for more info.

if exist('mipsver', 'file') == 2
    v = mipsver('all');
    vstr = v.Version;
    rdate = v.Date;
    TorF = 1;
else
    vstr = '';
    rdate = '';
    TorF = 0;
end
