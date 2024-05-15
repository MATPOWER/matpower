function rv = mostver(varargin)
% mostver - Prints or returns installed MOST version info.
% ::
%
%   mostver
%   v = mostver
%   v = mostver('all')
%
% When called with an output argument and no input argument, mostver
% returns the current MOST version numbers. With an input argument
% (e.g. ``'all'``) it returns  a struct with the fields ``Name``,
% ``Version``, ``Release``, and ``Date`` *(all char arrays)*. Calling
% mostver without assigning the return value prints the version and
% release date of the current installation of MOST.
%
% See also mpver.

%   MOST
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

v = struct( 'Name',     'MOST', ...
            'Version',  '1.3', ...
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
