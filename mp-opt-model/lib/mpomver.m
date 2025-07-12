function rv = mpomver(varargin)
% mpomver - Prints or returns installed MP-Opt-Model version info.
% ::
%
%   mpomver
%   v = mpomver
%   v = mpomver('all')
%
% When called with an output argument and no input argument, mpomver
% returns the current MP-Opt-Model version numbers. With an input argument
% (e.g. ``'all'``) it returns  a struct with the fields ``Name``,
% ``Version``, ``Release``, and ``Date`` *(all char arrays)*. Calling
% mpomver without assigning the return value prints the version and
% release date of the current installation of MP-Opt-Model.
%
% See also mpver.

%   MP-Opt-Model
%   Copyright (c) 2010-2025, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

v = struct( 'Name',     'MP-Opt-Model', ...
            'Version',  '5.0', ...
            'Release',  '', ...
            'Date',     '12-Jul-2025' );
if nargout > 0
    if nargin > 0
        rv = v;
    else
        rv = v.Version;
    end
else
    fprintf('%-22s Version %-9s  %11s\n', v.Name, v.Version, v.Date);
end
