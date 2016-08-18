function ef = cpf_nose_event(cb_data, cc)
%CPF_NOSE_EVENT  Event function to detect the nose point
%
%   EF = CPF_NOSE_EVENT(CB_DATA, CC)
%   
%   Inputs:
%       CB_DATA : struct of data for callback functions
%       CC : struct containing info about current point (continuation soln)
%
%   Outputs:
%       EF : event function value

%   MATPOWER
%   Copyright (c) 2016 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% event function value is dlam, the last element of the
%% normalized tangent vector at the current soln
ef = cc.z(end);

% fprintf('==== NOSE event fcn = %g\n', ef);
