function [TorF, vstr, rdate] = have_feature_glpk()
% have_feature_glpk - Detect availability/version info for GLPK.
%
% Private feature detection function implementing ``'glpk'`` tag for
% have_feature to detect availability/version of GLPK (GNU Linear
% Programming Kit), LP/MILP solver.
%
% See also have_feature, qps_master, miqps_master, glpk.

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
if exist('glpk','file') == 3    %% Windows OPTI install (no glpk.m)
    TorF = 1;
    str = evalc('glpk');
    pat = 'GLPK: GNU Linear Programming Kit \[v([^\s,\]]+).*\]';  %% OPTI, Giorgetti/Currie
    [s,e,tE,m,t] = regexp(str, pat);
    if ~isempty(t)
        vstr = t{1}{1};
    end
    pat = 'Built ([^\],])+';  %% OPTI, Giorgetti/Currie
    [s,e,tE,m,t] = regexp(str, pat);
    if ~isempty(t)
        rdate = datestr(t{1}{1}, 'dd-mmm-yyyy');
    end
elseif exist('glpk','file') == 2    %% others have glpk.m and ...
    if exist('__glpk__','file') == 3    %% octave __glpk__ MEX
        TorF = 1;
        if have_feature('evalc')
            str = evalc('glpk(1, 1, 1, 1, 1, ''U'', ''C'', -1, struct(''msglev'', 3))');
            pat = 'GLPK Simplex Optimizer, v([^\s,]+)';
            [s,e,tE,m,t] = regexp(str, pat);
            if ~isempty(t)
                vstr = t{1}{1};
            end
        end
    elseif exist('glpkcc','file') == 3  %% MATLAB glpkcc MEX
        TorF = 1;
        str = evalc('glpk');
        pat = 'GLPK Matlab interface\. Version: ([^\s,]+)';     %% glpkccm, Giorgetti/Klitgord
        [s,e,tE,m,t] = regexp(str, pat);
        if ~isempty(t)
            vstr = t{1}{1};
        end
    end
end
