function [TorF, vstr, rdate] = have_feature_rithmaticker()
% have_feature_rithmaticker - Detect availability/version info for rithmaticker.
%
% Private feature detection function implementing ``'rithmaticker'`` tag for
% have_feature to detect availability of rithmaticker (have_feature test
% function).
%
% See also have_feature.

%   MP-Test
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

TorF = exist('rithmaticker', 'file') == 2;
if TorF
    vstr = '3.1.4';
    rdate = datestr([2019 5 30 0 0 0], 'dd-mmm-yyyy');
else
    vstr = '';
    rdate = '';
end
