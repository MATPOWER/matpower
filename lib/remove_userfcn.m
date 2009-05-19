function mpc = remove_userfcn(mpc, stage, fcn)
%REMOVE_USERFCN Removes a userfcn from the list to be called for a case.
%
%   mpc = remove_userfcn(mpc, stage, fcn)
%
%   A userfcn is a callback function that can be called automatically by
%   MATPOWER at one of various stages in a simulation. This function removes
%   the last instance of the userfcn for the given stage with the name
%   specified by fcn.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

found = 0;
n = length(mpc.userfcn.(stage));
for k = n:-1:1
    if strcmp(mpc.userfcn.(stage)(k).fcn, fcn)
        found = 1;
        mpc.userfcn.(stage)(k) = [];
        break;
    end
end
if ~found
    error('remove_userfcn: userfcn ''%s'' for stage ''%s'' not found', fcn, stage);
end

return;
