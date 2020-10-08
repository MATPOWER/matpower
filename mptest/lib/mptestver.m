function rv = mptestver(varargin)
%MPTESTVER  Prints or returns MP-Test version info.
%   V = MPTESTVER returns the current MP-Test version numbers.
%   V = MPTESTVER('all') returns a struct with the fields Name, Version,
%   Release and Date (all char arrays). Calling MPTESTVER without assigning
%   the return value prints the version and release date of the current
%   installation of MP-Test.

%   MP-Test
%   Copyright (c) 2010-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

v = struct( 'Name',     'MP-Test', ... 
            'Version',  '7.1', ...
            'Release',  '', ...
            'Date',     '08-Oct-2020' );
if nargout > 0
    if nargin > 0
        rv = v;
    else
        rv = v.Version;
    end
else
    fprintf('%-22s Version %-10s %-11s\n', v.Name, v.Version, v.Date);
end
