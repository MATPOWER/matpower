function rv = gurobiver(varargin)
%GUROBIVER  Prints or returns GUROBI version info.
%   V = GUROBIVER returns the current GUROBI version numbers.
%   V = GUROBIVER('all') returns a struct with the fields Name, Version,
%   Release and Date (all strings). Calling GUROBIVER without assigning the
%   return value prints the version and release date of the current
%   installation of GUROBI.
%
%   See also MPVER.

%   MATPOWER
%   Copyright (c) 2010-2015 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   $Id$
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://matpower.org/ for more info.

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
    vn = sprintf('%d.%d.%d', result.versioninfo.major, result.versioninfo.minor, result.versioninfo.technical);
catch gurobiError
    fprintf('Gurobi Error!\n');
    disp(gurobiError.message);
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
    fprintf('%-22s Version %-10s %-11s\n', v.Name, v.Version, v.Date);
end
