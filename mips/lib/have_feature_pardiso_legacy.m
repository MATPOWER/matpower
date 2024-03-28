function [TorF, vstr, rdate] = have_feature_pardiso_legacy()
% have_feature_pardiso_legacy - Detect availability/version info for PARDISO (legacy interface).
%
% Private feature detection function implementing ``'pardiso_legacy'`` tag
% for have_feature to detect support for the legacy (v5.x) PARDISO interface,
% with individual MEX files for factor, solve, etc.
%
% See also have_feature, have_feature_pardiso, have_feature_pardiso_object.

%   MIPS
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MIPS.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mips for more info.

TorF = exist('pardisoinit', 'file') == 3 && ...
        exist('pardisoreorder', 'file') == 3 && ...
        exist('pardisofactor', 'file') == 3 && ...
        exist('pardisosolve', 'file') == 3 && ...
        exist('pardisofree', 'file') == 3;
rdate = '';
if TorF
    vstr = '6.x+';
    try
        A = sparse([1 2; 3 4]);
        b = [1;1];
        info = pardisoinit(11, 0);
        info = pardisoreorder(A, info, false);
%         % Summary PARDISO 5.1.0: ( reorder to reorder )
%         pat = 'Summary PARDISO (\.*\d)+:';
%         [s,e,tE,m,t] = regexp(evalc('info = pardisoreorder(A, info, true);'), pat);
%         if ~isempty(t)
%             vstr = t{1}{1};
%         end
        info = pardisofactor(A, info, false);
        [x, info] = pardisosolve(A, b, info, false);
        pardisofree(info);
        if any(x ~= [-1; 1])
            TorF = 0;
        end
    catch
        TorF = 0;
    end
else
    vstr = '';
end
