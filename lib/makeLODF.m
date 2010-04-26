function LODF = makeLODF(branch, PTDF);
%MAKELODF   Builds the line outage distribution factor matrix.
%   LODF = MAKELODF(BRANCH, PTDF) returns the DC line outage
%   distribution factor matrix for a given PTDF. The matrix is nbr x nbr,
%   where nbr is the number of branches.
%
%   Example:
%       H = makePTDF(baseMVA, bus, branch);
%       LODF = makeLODF(branch, H);
%
%   See also MAKEPTDF.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008-2010 by Power System Engineering Research Center (PSERC)
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

[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

[nl, nb] = size(PTDF);
f = branch(:, F_BUS);
t = branch(:, T_BUS);
Cft =  sparse([f; t], [1:nl 1:nl]', [ones(nl, 1); -ones(nl, 1)], nb, nl);

H = PTDF * Cft;
h = diag(H, 0);
LODF = H ./ (ones(nl, nl) - ones(nl, 1) * h');
LODF = LODF - diag(diag(LODF)) - eye(nl, nl);
