function Va = dcpf(B, Pbus, Va0, ref, pv, pq)
%DCPF  Solves a DC power flow.
%   [Va, success] = dcpf(B, Pbus, Va0, ref, pv, pq) solves for the bus
%   voltage angles at all but the reference bus, given the full system
%   B matrix and the vector of bus real power injections, the initial
%   vector of bus voltage angles (in radians), and column vectors with
%   the lists of bus indices for the swing bus, PV buses, and PQ buses,
%   respectively. Returns a vector of bus voltage angles in radians.

%   MATPOWER
%   $Id$
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   and Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% initialize result vector
Va = Va0;

%% update angles for non-reference buses
Va([pv; pq]) = B([pv; pq], [pv; pq]) \ (Pbus([pv; pq]) - B([pv; pq], ref) * Va0(ref));

return;
