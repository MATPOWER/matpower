function [TorF, vstr, rdate] = have_feature_evalc()
% have_feature_evalc - Detect availability/version info for :func:`evalc`.
%
% Private feature detection function implementing ``'evalc'`` tag for
% have_feature to detect support for :func:`evalc` function.
%
% See also have_feature, evalc.

%   MP-Opt-Model
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

TorF = exist('evalc') ~= 0;
if TorF
    if have_feature('matlab')
        v = have_feature('matlab', 'all');
    else
        v = have_feature('octave', 'all');
        TorF = 0;   %% In Octave, evalc was introduced in 4.2.x, but as of 6.x
                    %% we still choose to ignore evalc in Octave since it does
                    %% not capture the output of all functions (e.g. .mex/.oct
                    %% functions such as ipopt, glpk), which is one of the
                    %% main reasons we'd like to use it. Not only does it not
                    %% capture the output, it also does not prevent it from
                    %% going to the console.
    end
    vstr  = v.vstr;
    rdate = v.date;
else
    vstr = '';
    rdate = '';
end
