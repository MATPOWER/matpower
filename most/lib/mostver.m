function rv = mostver(varargin)
%MOSTVER  Prints or returns MOST version info for current installation.
%   V = MOSTVER returns the current MOST version number.
%   V = MOSTVER('all') returns a struct with the fields Name, Version,
%   Release and Date (all strings). Calling MOSTVER without assigning the
%   return value prints the version and release date of the current
%   installation of MOST.
%
%   See also MPVER.

%   MOST
%   Copyright (c) 2010-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

v = struct( 'Name',     'MOST', ... 
            'Version',  '1.0', ...
            'Release',  '', ...
            'Date',     '16-Dec-2016' );
if nargout > 0
    if nargin > 0
        rv = v;
    else
        rv = v.Version;
    end
else
    fprintf('%-22s Version %-9s  %11s\n', v.Name, v.Version, v.Date);
end
