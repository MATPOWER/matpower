function rv = mipsver(varargin)
%MIPSVER  Prints or returns MIPS version info for current installation.
%   V = MIPSVER returns the current MIPS version number.
%   V = MIPSVER('all') returns a struct with the fields Name, Version,
%   Release and Date (all strings). Calling MIPSVER without assigning the
%   return value prints the version and release date of the current
%   installation of MIPS.
%
%   See also MPVER.

%   MIPS
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2010 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MIPS.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MIPS is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MIPS is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MIPS. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MIPS, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MIPS grant
%   you additional permission to convey the resulting work.

v = struct( 'Name',     'MIPS', ... 
            'Version',  '1.0', ...
            'Release',  '', ...
            'Date',     '07-Feb-2011' );
if nargout > 0
    if nargin > 0
        rv = v;
    else
        rv = v.Version;
    end
else
    fprintf('%-22s Version %-9s  %11s\n', v.Name, v.Version, v.Date);
end
