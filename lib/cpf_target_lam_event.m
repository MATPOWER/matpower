function ef = cpf_target_lam_event(cb_data, cx)
%CPF_TARGET_LAM_EVENT  Event function to detect a target lambda value
%   EF = CPF_TARGET_LAM_EVENT(CB_DATA, CX)
%
%   CPF event function to detect the completion of the continuation curve
%   or another target value of lambda.
%
%   Inputs:
%       CB_DATA : struct of data for callback functions
%       CX : struct containing info about current point (continuation soln)
%
%   Outputs:
%       EF : event function value

%   MATPOWER
%   Copyright (c) 2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% event function value is a scalar equal to: current lambda - target lambda
target = cb_data.mpopt.cpf.stop_at;
if ischar(target)
    target = 0;
end
ef = cx.lam - target;
