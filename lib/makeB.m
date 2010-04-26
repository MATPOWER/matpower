function [Bp, Bpp] = makeB(baseMVA, bus, branch, alg)
%MAKEB   Builds the FDPF matrices, B prime and B double prime.
%   [BP, BPP] = MAKEB(BASEMVA, BUS, BRANCH, ALG) returns the two
%   matrices B prime and B double prime used in the fast decoupled power
%   flow. Does appropriate conversions to p.u. ALG is the value of the
%   PF_ALG option specifying the power flow algorithm.
%
%   Example:
%       [Bp, Bpp] = makeB(baseMVA, bus, branch, alg);
%
%   See also FDPF.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2010 by Power System Engineering Research Center (PSERC)
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

%% constants
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of lines

%% define named indices into bus, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%%-----  form Bp (B prime)  -----
temp_branch = branch;                       %% modify a copy of branch
temp_bus = bus;                             %% modify a copy of bus
temp_bus(:, BS) = zeros(nb, 1);             %% zero out shunts at buses
temp_branch(:, BR_B) = zeros(nl, 1);        %% zero out line charging shunts
temp_branch(:, TAP) = ones(nl, 1);          %% cancel out taps
if alg == 2                                 %% if XB method
    temp_branch(:, BR_R) = zeros(nl, 1);        %% zero out line resistance
end
Bp = -imag( makeYbus(baseMVA, temp_bus, temp_branch) );

%%-----  form Bpp (B double prime)  -----
if nargout == 2
    temp_branch = branch;                       %% modify a copy of branch
    temp_branch(:, SHIFT) = zeros(nl, 1);       %% zero out phase shifters
    if alg == 3                                 %% if BX method
        temp_branch(:, BR_R) = zeros(nl, 1);        %% zero out line resistance
    end
    Bpp = -imag( makeYbus(baseMVA, bus, temp_branch) );
end
