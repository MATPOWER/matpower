function mpc = get_mpc(om)
%GET_MPC  Returns the MATPOWER case struct.
%   MPC = GET_MPC(OM)
%
%   See also OPF_MODEL.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

mpc = om.mpc;
