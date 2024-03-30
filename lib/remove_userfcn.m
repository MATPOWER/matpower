function mpc = remove_userfcn(mpc, stage, fcn)
% remove_userfcn - Removes a userfcn from the list to be called for a case.
% ::
%
%   MPC = REMOVE_USERFCN(MPC, STAGE, FCN)
%
%   A userfcn is a callback function that can be called automatically by
%   MATPOWER at one of various stages in a simulation. This function removes
%   the last instance of the userfcn for the given STAGE with the function
%   handle specified by FCN.
%
% See also add_userfcn, run_userfcn, toggle_reserves, toggle_iflims,
% runopf_w_res.

%   MATPOWER
%   Copyright (c) 2009-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

n = length(mpc.userfcn.(stage));

if have_feature('octave')
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