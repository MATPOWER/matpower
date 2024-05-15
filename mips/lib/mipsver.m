function rv = mipsver(varargin)
% mipsver - Prints or returns installed MIPS version info.
% ::
%
%   mipsver
%   v = mipsver
%   v = mipsver('all')
%
% When called with an output argument and no input argument, mipsver
% returns the current MIPS version numbers. With an input argument
% (e.g. ``'all'``) it returns  a struct with the fields ``Name``,
% ``Version``, ``Release``, and ``Date`` *(all char arrays)*. Calling
% mipsver without assigning the return value prints the version and
% release date of the current installation of MIPS.
%
% See also mpver.

%   MIPS
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MIPS.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mips for more info.

v = struct( 'Name',     'MIPS', ... 
            'Version',  '1.5.1', ...
            'Release',  '', ...
            'Date',     '10-May-2024' );
if nargout > 0
    if nargin > 0
        rv = v;
    else
        rv = v.Version;
    end
else
    fprintf('%-22s Version %-9s  %11s\n', v.Name, v.Version, v.Date);
end
