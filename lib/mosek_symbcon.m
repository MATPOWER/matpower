function sc = mosek_symbcon
%MOSEK_SYMBCON  Returns struct containing MOSEK's symbolic constants
%   SC = MOSEK_SYMBCON
%
%   Returns a stuct containing all of MOSEk's symbolic constants, such
%   as those used to select the optimizer algorithm. Since the values
%   of these constants are not necessarily the same from one version of
%   MOSEK to the next, it is safer to use the symbolic constant rather
%   than the value in your MATPOWER code. For example, the following code
%
%       mpopt = mpoption('opf.dc.solver', 'MOSEK', 'mosek_lp_alg', 4);
%
%   would select primal simplex solver in MOSEK v6.x, but dual-simplex in
%   MOSEK v7.x. The recommended way to select a dual-simplex solver that
%   should work regardless of MOSEK version is to use the following:
%
%       sc = mosek_symbcon;
%       mpopt = mpoption('opf.dc.solver', 'MOSEK', ...
%                           'mosek_lp_alg', sc.MSK_OPTIMIZER_DUAL_SIMPLEX);

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2015 by Power System Engineering Research Center (PSERC)
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

[r, res] = mosekopt('symbcon echo(0)');
sc = res.symbcon;
