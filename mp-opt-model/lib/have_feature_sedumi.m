function [TorF, vstr, rdate] = have_feature_sedumi()
% have_feature_sedumi - Detect availability/version info for SeDuMi.
%
% Private feature detection function implementing ``'sedumi'`` tag for
% have_feature to detect availability/version of SeDuMi SDP solver
% (http://sedumi.ie.lehigh.edu).
%
% See also have_feature, sedumi.

%   MP-Opt-Model
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

TorF = exist('sedumi','file') == 2;
vstr = '';
rdate = '';
if TorF && have_feature('evalc')
    warn_state = warning;  %% sedumi turns (and leaves!) off all warnings
    str = evalc('x = sedumi([1 1], 1, [1;2])');
    warning(warn_state);
    pat = 'SeDuMi\s+([^\s]+)';
    [s,e,tE,m,t] = regexp(str, pat);
    if ~isempty(t)
        vstr = t{1}{1};
    end
end
