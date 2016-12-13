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
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010.
%             http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf

%   MATPOWER
%   Copyright (c) 2008-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define
nl = length(lam);

diaglam = sparse(1:nl, 1:nl, lam, nl, nl);
diagSbr_conj = sparse(1:nl, 1:nl, conj(Sbr), nl, nl);

[Saa, Sav, Sva, Svv] = d2Sbr_dV2(Cbr, Ybr, V, diagSbr_conj * lam);
Haa = 2 * real( Saa + dSbr_dVa.' * diaglam * conj(dSbr_dVa) );
Hva = 2 * real( Sva + dSbr_dVm.' * diaglam * conj(dSbr_dVa) );
Hav = 2 * real( Sav + dSbr_dVa.' * diaglam * conj(dSbr_dVm) );
Hvv = 2 * real( Svv + dSbr_dVm.' * diaglam * conj(dSbr_dVm) );
