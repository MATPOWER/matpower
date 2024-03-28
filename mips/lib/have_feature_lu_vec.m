function [TorF, vstr, rdate] = have_feature_lu_vec()
% have_feature_lu_vec - Detect availability/version info for LU vector support.
%
% Private feature detection function implementing ``'lu_vec'`` tag for
% have_feature to detect support for ``lu(..., 'vector')`` syntax.
%
%   See also have_feature, lu.

%   MIPS
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MIPS.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mips for more info.

vstr = '';
rdate = '';
if have_feature('matlab') && have_feature('matlab', 'vnum') < 7.003
    TorF = 0;     %% lu(..., 'vector') syntax not supported
else
    TorF = 1;
end
