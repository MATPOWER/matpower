function [TorF, vstr, rdate] = have_feature_octave()
%HAVE_FEATURE_OCTAVE  Detect availability/version info for Octave
%
%   Feature detection function implementing 'octave' tag for HAVE_FEATURE
%   to detect whether code is running under GNU Octave.
%
%   See also HAVE_FEATURE.

%   MP-Test
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

TorF = exist('OCTAVE_VERSION', 'builtin') == 5;
if TorF
    % v = ver('octave');
    %% workaround for https://savannah.gnu.org/bugs/index.php?59125
    v = ver;
    for k = 1:length(v)
        if strcmp(v(k).Name, 'Octave')
            v = v(k);
            break;
        end
    end

    vstr = v.Version;
    if ~isempty(v.Date)
        rdate = datestr(v.Date, 'dd-mmm-yyyy');
    else
        rdate = '';
    end
else
    vstr = '';
    rdate = '';
end
