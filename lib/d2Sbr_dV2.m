function [Haa, Hav, Hva, Hvv] = d2Sbr_dV2(Cbr, Ybr, V, lam)
%D2SBR_DV2   Computes 2nd derivatives of complex power flow w.r.t. voltage.
%   [HAA, HAV, HVA, HVV] = D2SBR_DV2(CBR, YBR, V, LAM) returns 4 matrices
%   containing the partial derivatives w.r.t. voltage angle and magnitude
%   of the product of a vector LAM with the 1st partial derivatives of the
%   complex branch power flows. Takes sparse connection matrix CBR, sparse
%   branch admittance matrix YBR, voltage vector V and nl x 1 vector of
%   multipliers LAM. Output matrices are sparse.
%
%   Example:
%       f = branch(:, F_BUS);
%       Cf =  sparse(1:nl, f, ones(nl, 1), nl, nb);
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       Cbr = Cf;
%       Ybr = Yf;
%       [Haa, Hav, Hva, Hvv] = d2Sbr_dV2(Cbr, Ybr, V, lam);
%
%   Here the output matrices correspond to:
%       Haa = (d/dVa (dSbr_dVa.')) * lam
%       Hav = (d/dVm (dSbr_dVa.')) * lam
%       Hva = (d/dVa (dSbr_dVm.')) * lam
%       Hvv = (d/dVm (dSbr_dVm.')) * lam

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
%   other modules (such as M-files and MEX-files) available in a
%   Matlab (or compatible) environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

%% define
nl = length(lam);
nb = length(V);

diaglam = sparse(1:nl, 1:nl, lam, nl, nl);
diagV   = sparse(1:nb, 1:nb, V, nb, nb);

A = Ybr' * diaglam * Cbr;
B = conj(diagV) * A * diagV;
D = sparse(1:nb, 1:nb, (A*V) .* conj(V), nb, nb);
E = sparse(1:nb, 1:nb, (A.'*conj(V)) .* V, nb, nb);
F = B + B.';
G = sparse(1:nb, 1:nb, ones(nb, 1)./abs(V), nb, nb);

Haa = F - D - E;
Hva = 1j * G * (B - B.' - D + E);
Hav = Hva.';
Hvv = G * F * G;
