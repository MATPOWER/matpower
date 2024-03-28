function rv = get_userdata(obj, name)
% get_userdata - Used to retrieve values of user data.
% ::
%
%   VAL = OBJ.GET_USERDATA(NAME) returns the value specified by the given name
%   or an empty matrix if userdata with NAME does not exist.
%
%   This function allows the user to retrieve any arbitrary data that was
%   saved in the object for later use. Data for a given NAME is saved by
%   assigning it to OBJ.userdata.(NAME).
%
%   This can be useful, for example, when using a user function to add
%   variables or constraints, etc. Suppose some special indexing is
%   constructed when adding some variables or constraints. This indexing data
%   can be stored and used later to "unpack" the results of the solved case.
%
% See also opt_model.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if isfield(obj.userdata, name)
    rv = obj.userdata.(name);
else
    rv = [];
end
