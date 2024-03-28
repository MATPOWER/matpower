function rv = gurobiver(varargin)
% gurobiver - Prints or returns installed GUROBI version info.
% ::
%
%   gurobiver
%   v = gurobiver
%   v = gurobiver('all')
%
% When called with an output argument and no input argument, gurobiver
% returns the current GUROBI version numbers. With an input argument
% (e.g. ``'all'``) it returns  a struct with the fields ``Name``,
% ``Version``, ``Release``, and ``Date`` *(all char arrays)*. Calling
% gurobiver without assigning the return value prints the version and
% release date of the current installation of GUROBI.
%
% See also mpver, gurobi.

%   MP-Opt-Model
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

g = have_feature('gurobi', 'all');
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
