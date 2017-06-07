function rv = get_userdata(om, name)
%GET_USERDATA  Used to retrieve values of user data.
%
%   VAL = OM.GET_USERDATA(NAME) returns the value specified by the given name
%
%   This function allows the user to retrieve any arbitrary data that was
%   saved in the object for later use. This can be useful when using a user
%   function to add variables, constraints, costs, etc. For example, suppose
%   some special indexing is constructed when adding some variables or
%   constraints. This indexing data can be stored and used later to "unpack"
%   the results of the solved case.
%
%   See also OPT_MODEL, SET_USERDATA.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if isfield(om.userdata, name)
    rv = om.userdata.(name);
else
    rv = [];
end
