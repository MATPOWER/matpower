function mpc = remove_userfcn(mpc, stage, fcn)
%REMOVE_USERFCN Removes a userfcn from the list to be called for a case.
%   MPC = REMOVE_USERFCN(MPC, STAGE, FCN)
%
%   A userfcn is a callback function that can be called automatically by
%   MATPOWER at one of various stages in a simulation. This function removes
%   the last instance of the userfcn for the given STAGE with the function
%   handle specified by FCN.
%
%   See also ADD_USERFCN, RUN_USERFCN, TOGGLE_RESERVES, TOGGLE_IFLIMS,
%   RUNOPF_W_RES.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

n = length(mpc.userfcn.(stage));

if have_fcn('octave')
    fcn_info = functions(fcn);
    for k = n:-1:1
        cb_info = functions(mpc.userfcn.(stage)(k).fcn);
        if strcmp(cb_info.function, fcn_info.function)
            mpc.userfcn.(stage)(k) = [];
            break;
        end
    end
else
    for k = n:-1:1
        if isequal(mpc.userfcn.(stage)(k).fcn, fcn)
            mpc.userfcn.(stage)(k) = [];
            break;
        end
    end
end