function mpc = get_mpc(om)
% get_mpc - Returns the MATPOWER case struct.
% ::
%
%   MPC = OM.GET_MPC()
%
% See also opt_model.

%   MATPOWER
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

mpc = om.mpc;
