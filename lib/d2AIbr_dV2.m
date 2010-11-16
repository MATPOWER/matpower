function [Haa, Hav, Hva, Hvv] = ...
    d2AIbr_dV2(dIbr_dVa, dIbr_dVm, Ibr, Ybr, V, lam)
%D2AIBR_DV2   Computes 2nd derivatives of |complex current|^2 w.r.t. V.
%   [HAA, HAV, HVA, HVV] = D2AIBR_DV2(DIBR_DVA, DIBR_DVM, IBR, YBR, V, LAM)
%   returns 4 matrices containing the partial derivatives w.r.t. voltage
%   angle and magnitude of the product of a vector LAM with the 1st partial
%   derivatives of the square of the magnitude of the branch currents.
%   Takes sparse first derivative matrices of complex flow, complex flow
%   vector, sparse branch admittance matrix YBR, voltage vector V and
%   nl x 1 vector of multipliers LAM. Output matrices are sparse.
%
%   Example:
%       f = branch(:, F_BUS);
%       Cf =  sparse(1:nl, f, ones(nl, 1), nl, nb);
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dIf_dVa, dIf_dVm, dIt_dVa, dIt_dVm, If, It] = ...
%               dIbr_dV(branch, Yf, Yt, V);
%       Cbr = Cf;
%       Ybr = Yf;
%       dIbr_dVa = dIf_dVa;
%       dIbr_dVm = dIf_dVm;
%       Ibr = If;
%       [Haa, Hav, Hva, Hvv] = ...
%             d2AIbr_dV2(dIbr_dVa, dIbr_dVm, Ibr, Ybr, V, lam);
%
%   Here the output matrices correspond to:
%     Haa = (d/dVa (dAIbr_dVa.')) * lam
%     Hav = (d/dVm (dAIbr_dVa.')) * lam
%     Hva = (d/dVa (dAIbr_dVm.')) * lam
%     Hvv = (d/dVm (dAIbr_dVm.')) * lam
%
%   See also DIBR_DV.
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
nl = length(lam);

diaglam = sparse(1:nl, 1:nl, lam, nl, nl);
diagIbr_conj = sparse(1:nl, 1:nl, conj(Ibr), nl, nl);

[Iaa, Iav, Iva, Ivv] = d2Ibr_dV2(Ybr, V, diagIbr_conj * lam);
Haa = 2 * real( Iaa + dIbr_dVa.' * diaglam * conj(dIbr_dVa) );
Hva = 2 * real( Iva + dIbr_dVm.' * diaglam * conj(dIbr_dVa) );
Hav = 2 * real( Iav + dIbr_dVa.' * diaglam * conj(dIbr_dVm) );
Hvv = 2 * real( Ivv + dIbr_dVm.' * diaglam * conj(dIbr_dVm) );
