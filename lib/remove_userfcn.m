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
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009-2010 by Power System Engineering Research Center (PSERC)
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