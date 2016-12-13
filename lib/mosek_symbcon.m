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
%   Copyright (c) 2015-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[r, res] = mosekopt('symbcon echo(0)');
sc = res.symbcon;
