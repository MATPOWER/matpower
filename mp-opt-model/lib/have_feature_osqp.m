function [TorF, vstr, rdate] = have_feature_osqp()
%HAVE_FEATURE_OSQP  Detect availability/version info for OSQP
%
%   Feature detection function implementing 'osqp' tag for HAVE_FEATURE
%   to detect availability/version of OSQP, (Operator Splitting QP solver)
%   (https://osqp.org)
%
%   See also HAVE_FEATURE, QPS_MASTER, OSQP.

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

TorF = exist('osqp', 'file') == 2;
vstr = '';
rdate = '';
if TorF
    try
        o = osqp();
        vstr = o.version();
    catch
        TorF = 0;
        fprintf('OSQP Error!\n');
    end
end
