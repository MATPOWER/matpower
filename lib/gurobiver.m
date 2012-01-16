function rv = gurobiver(varargin)
%GUROBIVER  Prints or returns GUROBI/GUROBI_MEX version info.
%   V = GUROBIVER returns the current GUROBI/GUROBI_MEX version numbers.
%   V = GUROBIVER('all') returns a struct with the fields Name, Version,
%   Release and Date (all strings). Calling GUROBIVER without assigning the
%   return value prints the version and release date of the current
%   installation of GUROBI.
%
%   See also MPVER.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2010, 2012 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MIPS.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MIPS is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MIPS is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MIPS. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MIPS, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MIPS grant
%   you additional permission to convey the resulting work.

[x, f, e, o] = gurobi_mex(1, 1, [], [], [], 1, 1, 'C', ...
	struct('Display', 0, 'DisplayInterval', Inf));
if isfield(o, 'Versions')
	vn = sprintf('%d.%d.%d/%.2f', o.Versions.LibMajor, ...
		o.Versions.LibMinor, o.Versions.LibTechi, o.Versions.GurobiMex);
elseif isfield(o, 'Lib.Ver.Major')
	vn = sprintf('%d.%d.%d/%.2f', o.('Lib.Ver.Major'), ...
		o.('Lib.Ver.Minor'), o.('Lib.Ver.Techi'), o.('GurobiMex.Ver'));
else
	vn = '<unknown>';
end

v = struct( 'Name',     'Gurobi', ... 
            'Version',  vn, ...
            'Release',  '', ...
            'Date',     '' );
if nargout > 0
    if nargin > 0
        rv = v;
    else
        rv = v.Version;
    end
else
    fprintf('%-22s Version %-10s %-11s   %s\n', v.Name, v.Version, v.Date, computer);
end
