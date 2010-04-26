function [Haa, Hav, Hva, Hvv] = ...
    d2ASbr_dV2(dSbr_dVa, dSbr_dVm, Sbr, Cbr, Ybr, V, lam)
%D2ASBR_DV2   Computes 2nd derivatives of |complex power flow|^2 w.r.t. V.
%   [HAA, HAV, HVA, HVV] = D2ASBR_DV2(DSBR_DVA, DSBR_DVM, SBR, CBR, YBR, V, LAM)
%   returns 4 matrices containing the partial derivatives w.r.t. voltage
%   angle and magnitude of the product of a vector LAM with the 1st partial
%   derivatives of the square of the magnitude of branch complex power flows.
%   Takes sparse first derivative matrices of complex flow, complex flow
%   vector, sparse connection matrix CBR, sparse branch admittance matrix YBR,
%   voltage vector V and nl x 1 vector of multipliers LAM. Output matrices
%   are sparse.
%
%   Example:
%       f = branch(:, F_BUS);
%       Cf =  sparse(1:nl, f, ones(nl, 1), nl, nb);
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = ...
%               dSbr_dV(branch, Yf, Yt, V);
%       Cbr = Cf;
%       Ybr = Yf;
%       dSbr_dVa = dSf_dVa;
%       dSbr_dVm = dSf_dVm;
%       Sbr = Sf;
%       [Haa, Hav, Hva, Hvv] = ...
%             d2ASbr_dV2(dSbr_dVa, dSbr_dVm, Sbr, Cbr, Ybr, V, lam);
%
%   Here the output matrices correspond to:
%     Haa = (d/dVa (dASbr_dVa.')) * lam
%     Hav = (d/dVm (dASbr_dVa.')) * lam
%     Hva = (d/dVa (dASbr_dVm.')) * lam
%     Hvv = (d/dVm (dASbr_dVm.')) * lam
%
%   See also DSBR_DV.

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
diagSbr_conj = sparse(1:nl, 1:nl, conj(Sbr), nl, nl);

[Saa, Sav, Sva, Svv] = d2Sbr_dV2(Cbr, Ybr, V, diagSbr_conj * lam);
Haa = 2 * real( Saa + dSbr_dVa.' * diaglam * conj(dSbr_dVa) );
Hva = 2 * real( Sva + dSbr_dVm.' * diaglam * conj(dSbr_dVa) );
Hav = 2 * real( Sav + dSbr_dVa.' * diaglam * conj(dSbr_dVm) );
Hvv = 2 * real( Svv + dSbr_dVm.' * diaglam * conj(dSbr_dVm) );
