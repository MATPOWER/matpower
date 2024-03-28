function [TorF, vstr, rdate] = have_feature_pardiso_object()
% have_feature_pardiso_object - Detect availability/version info for PARDISO (object interface).
%
% Private feature detection function implementing ``'pardiso_object'`` tag for
% have_feature to detect support for the object-oriented (v6.x and later)
% PARDISO interface.
%
% See also have_feature, have_feature_pardiso, have_feature_pardiso_legacy.

%   MIPS
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MIPS.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mips for more info.

TorF = exist('pardiso', 'file') == 2;
rdate = '';
if TorF
    vstr = '6.x+';
    try
        id = 1;
        A = sparse([1 2; 3 4]);
        b = [1;1];
        p = pardiso(id, 11, 0);
        p.factorize(id, A);
        x = p.solve(id, A, b);
        p.free(id);
        p.clear();
        if any(x ~= [-1; 1])
            TorF = 0;
        end
    catch
        TorF = 0;
    end
else
    vstr = '';
end
