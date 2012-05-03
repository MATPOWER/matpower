function rv = userdata(om, name, val)
%USERDATA  Used to save or retrieve values of user data.
%
%   OM = USERDATA(OM, NAME, VAL) saves the value under the given name.
%   VAL = USERDATA(OM, NAME) returns the value specified by the given name
%
%   This function allows the user to save any arbitrary data in the object
%   for later use. This can be useful when using a user function to add
%   variables, constraints, costs, etc. For example, suppose some special
%   indexing is constructed when adding some variables or constraints.
%   This indexing data can be stored and used later to "unpack" the results
%   of the solved case.
%
%   See also OPT_MODEL.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008-2012 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

if nargin == 3
    om.userdata.(name) = val;
    rv = om;
else
    if isfield(om.userdata, name)
        rv = om.userdata.(name);
    else
        rv = [];
    end
end
