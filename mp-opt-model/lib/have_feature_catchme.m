function [TorF, vstr, rdate] = have_feature_catchme()
%HAVE_FEATURE_CATCHME  Detect availability/version info for CATCH(ME) syntax.
%
%   Feature detection function implementing 'catchme' tag for HAVE_FEATURE
%   to detect support for CATCH ME syntax in TRY/CATCH constructs.
%
%   See also HAVE_FEATURE, TRY, CATCH.

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

vstr = '';
rdate = '';
if have_feature('octave')
    v = have_feature('octave', 'all');
    if v.vnum > 3.006
        TorF = 1;
        vstr = v.vstr;
        rdate = v.date;
    end
else
    v = have_feature('matlab', 'all');
    if v.vnum > 7.004
        TorF = 1;
        vstr = v.vstr;
        rdate = v.date;
    end
end
