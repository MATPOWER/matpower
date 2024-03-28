function [TorF, vstr, rdate] = have_feature_cplex()
% have_feature_cplex - Detect availability/version info for CPLEX.
%
% Private feature detection function implementing ``'cplex'`` tag for
% have_feature to detect availability/version of CPLEX (IBM ILOG CPLEX
% Optimizer).
%
% See also have_feature, qps_master, miqps_master, cplex, cplexqp, cplexlp.

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
if exist('cplexqp', 'file')
    %% it's installed, but we need to check for MEX for this arch
    p = which('cplexqp');   %% get the path
    len = length(p) - length('cplexqp.p');
    w = what(p(1:len));             %% look for mex files on the path
    for k = 1:length(w.mex)
        if regexp(w.mex{k}, 'cplexlink[^\.]*');
            TorF = 1;
            break;
        end
    end
end
if TorF
    try
        cplex = Cplex('null');
        vstr = cplex.getVersion;
    catch
        TorF = 0;
    end
end
