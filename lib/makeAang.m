function [Aang, lang, uang, iang]  = makeAang(baseMVA, branch, nb, mpopt)
%MAKEAANG  Construct constraints for branch angle difference limits.
%   [AANG, LANG, UANG, IANG]  = MAKEAANG(BASEMVA, BRANCH, NB, MPOPT)
%
%   Constructs the parameters for the following linear constraint limiting
%   the voltage angle differences across branches, where Va is the vector
%   of bus voltage angles. NB is the number of buses.
%
%       LANG <= AANG * Va <= UANG
%
%   IANG is the vector of indices of branches with angle difference limits.
%   The limits are given in the ANGMIN and ANGMAX columns of the branch
%   matrix. Voltage angle differences are taken to be unbounded below if
%   ANGMIN < -360 and unbounded above if ANGMAX > 360. If both ANGMIN and
%   ANGMAX are zero, the angle difference is assumed to be unconstrained.
%
%   Example:
%       [Aang, lang, uang, iang]  = makeAang(baseMVA, branch, nb, mpopt);

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
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

%% options
ignore_ang_lim = mpopt(25);     %% OPF_IGNORE_ANG_LIM

%% define named indices into data matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

if ignore_ang_lim
  Aang  = sparse(0, nb);
  lang  = [];
  uang  = [];
  iang  = [];
else
  iang = find((branch(:, ANGMIN) & branch(:, ANGMIN) > -360) | ...
              (branch(:, ANGMAX) & branch(:, ANGMAX) < 360));
  nang = length(iang);

  if nang > 0
    ii = [(1:nang)'; (1:nang)'];
    jj = [branch(iang, F_BUS); branch(iang, T_BUS)];
    Aang = sparse(ii, jj, [ones(nang, 1); -ones(nang, 1)], nang, nb);
    lang = branch(iang, ANGMIN) * pi/180;
    uang = branch(iang, ANGMAX) * pi/180;
  else
    Aang = sparse(0, nb);
    lang =[];
    uang =[];
  end
end
