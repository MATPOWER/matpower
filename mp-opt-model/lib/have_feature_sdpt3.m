function [TorF, vstr, rdate] = have_feature_sdpt3()
%HAVE_FEATURE_SDPT3  Detect availability/version info for SDPT3
%
%   Feature detection function implementing 'sdpt3' tag for HAVE_FEATURE
%   to detect availability/version of SDPT3 SDP solver
%   (https://github.com/sqlp/sdpt3)
%
%   See also HAVE_FEATURE, SDPT3.

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

TorF = exist('sdpt3','file') == 2;
vstr = '';
rdate = '';
if TorF && have_feature('evalc')
    str = evalc('help sdpt3');
    pat = 'version\s+([^\s]+).*Last Modified: ([^\n]+)\n';
    [s,e,tE,m,t] = regexp(str, pat);
    if ~isempty(t)
        vstr = t{1}{1};
        rdate = datestr(t{1}{2}, 'dd-mmm-yyyy');
    end
end
