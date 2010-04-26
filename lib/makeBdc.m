function [Bbus, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch)
%MAKEBDC   Builds the B matrices and phase shift injections for DC power flow.
%   [BBUS, BF, PBUSINJ, PFINJ] = MAKEBDC(BASEMVA, BUS, BRANCH) returns the
%   B matrices and phase shift injection vectors needed for a DC power flow.
%   The bus real power injections are related to bus voltage angles by
%       P = BBUS * Va + PBUSINJ
%   The real power flows at the from end the lines are related to the bus
%   voltage angles by
%       Pf = BF * Va + PFINJ
%   Does appropriate conversions to p.u.
%
%   Example:
%       [Bbus, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);
%
%   See also DCPF.

%   MATPOWER
%   $Id$
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   and Ray Zimmerman, PSERC Cornell
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

%% check that bus numbers are equal to indices to bus (one set of bus numbers)
if any(bus(:, BUS_I) ~= (1:nb)')
    error('makeBdc: buses must be numbered consecutively in bus matrix')
end

%% for each branch, compute the elements of the branch B matrix and the phase
%% shift "quiescent" injections, where
%%
%%      | Pf |   | Bff  Bft |   | Vaf |   | Pfinj |
%%      |    | = |          | * |     | + |       |
%%      | Pt |   | Btf  Btt |   | Vat |   | Ptinj |
%%
stat = branch(:, BR_STATUS);                    %% ones at in-service branches
b = stat ./ branch(:, BR_X);                    %% series susceptance
tap = ones(nl, 1);                              %% default tap ratio = 1
i = find(branch(:, TAP));                       %% indices of non-zero tap ratios
tap(i) = branch(i, TAP);                        %% assign non-zero tap ratios
b = b ./ tap;

%% build connection matrix Cft = Cf - Ct for line and from - to buses
f = branch(:, F_BUS);                           %% list of "from" buses
t = branch(:, T_BUS);                           %% list of "to" buses
i = [(1:nl)'; (1:nl)'];                         %% double set of row indices
Cft = sparse(i, [f;t], [ones(nl, 1); -ones(nl, 1)], nl, nb);    %% connection matrix

%% build Bf such that Bf * Va is the vector of real branch powers injected
%% at each branch's "from" bus
Bf = sparse(i, [f; t], [b; -b]);    % = spdiags(b, 0, nl, nl) * Cft;

%% build Bbus
Bbus = Cft' * Bf;

%% build phase shift injection vectors
Pfinj = b .* (-branch(:, SHIFT) * pi/180);      %% injected at the from bus ...
    % Ptinj = -Pfinj;                           %% ... and extracted at the to bus
Pbusinj = Cft' * Pfinj;                         %% Pbusinj = Cf * Pfinj + Ct * Ptinj;
