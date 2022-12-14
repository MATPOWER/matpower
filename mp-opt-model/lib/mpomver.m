function rv = mpomver(varargin)
%MPOMVER  Prints or returns MP-Opt-Model version info for current installation.
%   V = MPOMVER returns the current MP-Opt-Model version number.
%   V = MPOMVER('all') returns a struct with the fields Name, Version,
%   Release and Date (all strings). Calling MPOMVER without assigning the
%   return value prints the version and release date of the current
%   installation of MP-Opt-Model.
%
%   See also MPVER.

%   MP-Opt-Model
%   Copyright (c) 2010-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

v = struct( 'Name',     'MP-Opt-Model', ... 
            'Version',  '4.1', ...
            'Release',  '', ...
            'Date',     '13-Dec-2022' );
if nargout > 0
    if nargin > 0
        rv = v;
    else
        rv = v.Version;
    end
else
    fprintf('%-22s Version %-9s  %11s\n', v.Name, v.Version, v.Date);
end
