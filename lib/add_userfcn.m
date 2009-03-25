function mpc = add_userfcn(mpc, name, args, allow_multiple)
%ADD_USERFCN   Appends a userfcn to the list to be called for a case.
%   mpc = add_userfcn(mpc, name, args)
%   mpc = add_userfcn(mpc, name, args, allow_multiple)
%
%   mpc  : the case struct
%   name : the name of the userfcn
%   args : the value to be passed as an argument to the userfcn
%          (typically a struct)
%   allow_multiple : (optional) if TRUE, allows the same function to
%          be added more than once.
%
%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 4
    allow_multiple = 0;
end
if isfield(mpc, 'userfcn')
    n = length(mpc.userfcn) + 1;
    if ~allow_multiple
        for k = 1:n-1
            if strcmp(mpc.userfcn(k).name, name)
                error('add_userfcn: the function ''%s'' has already been added', name);
            end
        end
    end
else
    n = 1;
end

mpc.userfcn(n).name = name;
mpc.userfcn(n).args = args;

return;
