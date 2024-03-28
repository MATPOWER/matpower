function [TorF, vstr, rdate] = have_feature_fmincon()
% have_feature_fmincon - Detect availability/version info for :func:`fmincon`.
%
% Private feature detection function implementing ``'fmincon'`` tag for
% have_feature to detect availability/version of :func:`fmincon` from the
% MATLAB Optimization Toolbox.
%
% See also have_feature, nlps_master, fmincon.

%   MP-Opt-Model
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

TorF = 0;
vstr = '';
rdate = '';
v = have_feature('optim', 'all');
if v.av && have_feature('matlab')   %% ignore Octave version
    TorF = (exist('fmincon', 'file') == 2 || exist('fmincon', 'file') == 6);
    if TorF
        vstr = v.vstr;
        rdate = v.date;
    end
end
