function [TorF, vstr, rdate] = have_feature_matlab()
% have_feature_matlab - Detect availability/version info for MATLAB.
%
% Private feature detection function implementing ``'matlab'`` tag for
% have_feature to detect whether code is running under MATLAB.
%
%   See also have_feature.

%   MP-Test
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

v = ver('matlab');
if length(v) > 1
    warning('The built-in VER command is behaving strangely, probably as a result of installing a 3rd party toolbox in a directory named ''matlab'' on your path. Check each element of the output of ver(''matlab'') to find the offending toolbox, then move the toolbox to a more appropriately named directory.');
    v = v(1);
end
if ~isempty(v) && isfield(v, 'Version') && ~isempty(v.Version)
    TorF = 1;
    vstr = v.Version;
    rdate = v.Date;
else
    TorF = 0;
    vstr = '';
    rdate = '';
end
