function rv = mptestver(varargin)
%mptestver - Prints or returns installed MP-Test version info.
% ::
%
%   mptestver
%   v = mptestver
%   v = mptestver('all')
%
% When called with an output argument and no input argument, mptestver
% returns the current MP-Test version numbers. With an input argument
% (e.g. ``'all'``) it returns  a struct with the fields ``Name``,
% ``Version``, ``Release``, and ``Date`` *(all char arrays)*. Calling
% mptestver without assigning the return value prints the version and
% release date of the current installation of MP-Test.

%   MP-Test
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

v = struct( 'Name',     'MP-Test', ... 
            'Version',  '8.0', ...
            'Release',  '', ...
            'Date',     '10-May-2024' );
if nargout > 0
    if nargin > 0
        rv = v;
    else
        rv = v.Version;
    end
else
    fprintf('%-22s Version %-10s %-11s\n', v.Name, v.Version, v.Date);
end
