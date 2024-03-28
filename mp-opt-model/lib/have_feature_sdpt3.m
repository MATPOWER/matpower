function [TorF, vstr, rdate] = have_feature_sdpt3()
% have_feature_sdpt3 -  Detect availability/version info for SDPT3.
%
% Private feature detection function implementing ``'sdpt3'`` tag for
% have_feature to detect availability/version of SDPT3 SDP solver
% (https://github.com/sqlp/sdpt3)
%
% See also have_feature, sdpt3.

%   MP-Opt-Model
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
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
