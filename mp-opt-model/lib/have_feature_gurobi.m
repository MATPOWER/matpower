function [TorF, vstr, rdate] = have_feature_gurobi()
%HAVE_FEATURE_GUROBI  Detect availability/version info for Gurobi
%
%   Feature detection function implementing 'gurobi' tag for HAVE_FEATURE
%   to detect availability/version of Gurobi optimizer
%   (https://www.gurobi.com).
%
%   See also HAVE_FEATURE, QPS_MASTER, MIQPS_MASTER, GUROBI.

%   MP-Opt-Model
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

TorF = exist('gurobi', 'file') == 3;
vstr = '';
rdate = '';
if TorF
    try
        model = struct( ...
            'A', sparse(1), ...
            'rhs', 1, ...
            'sense', '=', ...
            'vtype', 'C', ...
            'obj', 1, ...
            'modelsense', 'min' ...
        );
        params = struct( ...
            'outputflag', 0 ...
        );
        result = gurobi(model, params);
        vstr = sprintf('%d.%d.%d', result.versioninfo.major, result.versioninfo.minor, result.versioninfo.technical);
    catch % gurobiError
        TorF = 0;
        fprintf('Gurobi Error!\n');
%         disp(gurobiError.message);
    end
end
