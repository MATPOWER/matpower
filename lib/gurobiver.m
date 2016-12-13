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
%   Copyright (c) 2010-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

g = have_fcn('gurobi', 'all');
if ~g.av
    g.vstr = '<unknown>';
end

v = struct( 'Name',     'Gurobi', ... 
            'Version',  g.vstr, ...
            'Release',  '', ...
            'Date',     g.date );
if nargout > 0
    if nargin > 0
        rv = v;
    else
        rv = v.Version;
    end
else
    if g.av
        fprintf('%-22s Version %-10s %-11s\n', v.Name, v.Version, v.Date);
    else
        fprintf('%-22s -- not installed --\n', v.Name);
    end
end
