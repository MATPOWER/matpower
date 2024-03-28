function [TorF, vstr, rdate] = have_feature_intlinprog()
% have_feature_intlinprog - Detect availability/version info for :func:`intlinprog`.
%
% Private feature detection function implementing ``'intlinprog'`` tag for
% have_feature to detect availability/version of :func:`intlinprog`, MILP
% solver from MATLAB Optimization Toolbox 7.0 (R2014a) and later.
%
% See also have_feature, miqps_master, intlinprog.

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
    TorF = exist('intlinprog', 'file') == 2;
    if TorF
        vstr = v.vstr;
        rdate = v.date;
    end
end
