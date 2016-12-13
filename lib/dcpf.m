function [Va, success] = dcpf(B, Pbus, Va0, ref, pv, pq)
%DCPF  Solves a DC power flow.
%   [VA, SUCCESS] = DCPF(B, PBUS, VA0, REF, PV, PQ) solves for the bus
%   voltage angles at all but the reference bus, given the full system
%   B matrix and the vector of bus real power injections, the initial
%   vector of bus voltage angles (in radians), and column vectors with
%   the lists of bus indices for the swing bus, PV buses, and PQ buses,
%   respectively. Returns a vector of bus voltage angles in radians.
%
%   See also RUNDCPF, RUNPF.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% constant
Va_threshold = 1e5;     %% arbitrary threshold on |Va| for declaring failure

%% initialize result vector
Va = Va0;
success = 1;    %% successful by default

%% set up to trap non-singular matrix warnings
[lastmsg, lastid] = lastwarn;
lastwarn('');

%% update angles for non-reference buses
Va([pv; pq]) = B([pv; pq], [pv; pq]) \ ...
                (Pbus([pv; pq]) - B([pv; pq], ref) * Va0(ref));

[msg, id] = lastwarn;
%% Octave is not consistent in assigning proper warning id, so we'll just
%% check for presence of *any* warning
if ~isempty(msg) || max(abs(Va)) > Va_threshold
    success = 0;
end

%% restore warning state
lastwarn(lastmsg, lastid);
