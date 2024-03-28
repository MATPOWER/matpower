function [TorF, vstr, rdate] = have_feature_linprog()
% have_feature_linprog - Detect availability/version info for :func:`linprog`.
%
% Private feature detection function implementing ``'linprog'`` tag for
% have_feature to detect availability/version of :func:`linprog`, LP
% solver from the MATLAB Optimization Toolbox.
%
% See also have_feature, have_feature_linprog_ds, qps_master, linprog.

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
    TorF = exist('linprog', 'file') == 2;
    if TorF
        vstr = v.vstr;
        rdate = v.date;
    end
end
