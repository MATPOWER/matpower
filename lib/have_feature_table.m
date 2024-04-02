function [TorF, vstr, rdate] = have_feature_table()
% have_feature_table - Detect availability/version info for :class:`table`.
%
% Private feature detection function implementing ``'table'`` tag for
% have_feature to detect availability/version of TABLE, included in
% MATLAB R2013b and as of this writing in Mar 2024, available for Octave
% as Tablicious: https://github.com/apjanke/octave-tablicious
%
% See also have_feature, table.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

TorF = exist('table', 'file') == 2;
vstr = '';
rdate = '';
if TorF
    if have_feature('matlab')
        vstr = have_feature('matlab', 'vstr');
        rdate = have_feature('matlab', 'date');
    else
        v = ver('tablicious');
        if ~isempty(v) && isfield(v, 'Version')
            vstr = v.Version;
        end
        if ~isempty(v) && isfield(v, 'Date')
            rdate = v.Date;
        end
    end
end
