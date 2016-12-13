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
%   Copyright (c) 2008-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

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
