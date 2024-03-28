function [TorF, vstr, rdate] = have_feature_gurobi()
% have_feature_gurobi - Detect availability/version info for GUROBI.
%
% Private feature detection function implementing ``'gurobi'`` tag for
% have_feature to detect availability/version of Gurobi optimizer
% (https://www.gurobi.com).
%
% See also have_feature, qps_master, miqps_master, gurobi.

%   MP-Opt-Model
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
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
