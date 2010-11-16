function [Haa, Hav, Hva, Hvv] = d2Ibr_dV2(Ybr, V, lam)
%D2IBR_DV2   Computes 2nd derivatives of complex branch current w.r.t. voltage.
%   [HAA, HAV, HVA, HVV] = D2IBR_DV2(CBR, YBR, V, LAM) returns 4 matrices
%   containing the partial derivatives w.r.t. voltage angle and magnitude
%   of the product of a vector LAM with the 1st partial derivatives of the
%   complex branch currents. Takes sparse branch admittance matrix YBR,
%   voltage vector V and nl x 1 vector of multipliers LAM. Output matrices
%   are sparse.
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       Ybr = Yf;
%       [Haa, Hav, Hva, Hvv] = d2Ibr_dV2(Ybr, V, lam);
%
%   Here the output matrices correspond to:
%       Haa = (d/dVa (dIbr_dVa.')) * lam
%       Hav = (d/dVm (dIbr_dVa.')) * lam
%       Hva = (d/dVa (dIbr_dVm.')) * lam
%       Hvv = (d/dVm (dIbr_dVm.')) * lam
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010.
%             http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf

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

%% define
nb = length(V);

diaginvVm = sparse(1:nb, 1:nb, ones(nb, 1)./abs(V), nb, nb);

Haa = sparse(1:nb, 1:nb, -(Ybr.' * lam) .* V, nb, nb);
Hva = -1j * Haa * diaginvVm;
Hav = Hva;
Hvv = sparse(nb, nb);
