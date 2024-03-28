function rv = osqpver(varargin)
% osqpver - Prints or returns installed OSQP version info.
% ::
%
%   osqpver
%   v = osqpver
%   v = osqpver('all')
%
% When called with an output argument and no input argument, osqpver
% returns the current OSQP version numbers. With an input argument
% (e.g. ``'all'``) it returns  a struct with the fields ``Name``,
% ``Version``, ``Release``, and ``Date`` *(all char arrays)*. Calling
% osqpver without assigning the return value prints the version and
% release date of the current installation of OSQP.
%
% See also mpver, osqp.

%   MP-Opt-Model
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

o = have_feature('osqp', 'all');
if ~o.av
    o.vstr = '<unknown>';
end

v = struct( 'Name',     'OSQP', ... 
            'Version',  o.vstr, ...
            'Release',  '', ...
            'Date',     o.date );
if nargout > 0
    if nargin > 0
        rv = v;
    else
        rv = v.Version;
    end
else
    if o.av
        fprintf('%-22s Version %-10s %-11s\n', v.Name, v.Version, v.Date);
    else
        fprintf('%-22s -- not installed --\n', v.Name);
    end
end
