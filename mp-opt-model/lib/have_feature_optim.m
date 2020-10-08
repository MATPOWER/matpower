function [TorF, vstr, rdate] = have_feature_optim()
%HAVE_FEATURE_OPTIM  Detect availability/version info for Optimization Toolbox
%
%   Feature detection function implementing 'optim' tag for HAVE_FEATURE
%   to detect availability/version of the Optimization Toolbox.
%
%   See also HAVE_FEATURE.

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
matlab = have_feature('matlab');
if ~matlab || (matlab && license('test', 'optimization_toolbox'))
    v = ver('optim');
    if length(v) > 1
        warning('The built-in VER command is behaving strangely, probably as a result of installing a 3rd party toolbox in a directory named ''optim'' on your path. Check each element of the output of ver(''optim'') to find the offending toolbox, then move the toolbox to a more appropriately named directory.');
        v = v(1);
    end
    if ~isempty(v) && ~isempty(v.Version)
        TorF = 1;
        vstr = v.Version;
        rdate = v.Date;
    end
end
