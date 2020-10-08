function [TorF, vstr, rdate] = have_feature_quadprog()
%HAVE_FEATURE_QUADPROG  Detect availability/version info for QUADPROG
%
%   Feature detection function implementing 'quadprog' tag for HAVE_FEATURE
%   to detect availability/version of QUADPROG, QP solver from the MATLAB
%   Optimization Toolbox.
%
%   See also HAVE_FEATURE, HAVE_FEATURE_QUADPROG_LS, QPS_MASTER, QUADPROG.

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
v = have_feature('optim', 'all');
if v.av
    TorF = exist('quadprog', 'file') == 2;
    %% Octave optim 1.5.0 and earlier, had problems with
    %% incorrect lambdas, including opposite sign
    %% convention for equality multipliers
    if have_feature('octave') && v.vnum <= 1.005
        TorF = 0;
    end
    if TorF
        vstr = v.vstr;
        rdate = v.date;
    end
end
