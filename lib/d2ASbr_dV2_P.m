function [Haa, Hav, Hva, Hvv] = ...
    d2ASbr_dV2_P(dSbr_dVa, dSbr_dVm, Sbr, Cbr, Ybr, V, mu)
%d2ASbr_dV2_P   Computes 2nd derivatives of |complex power flow|^2 w.r.t. V.
%   [HAA, HAV, HVA, HVV] = d2ASbr_dV2_P(DSBR_DVA, DSBR_DVM, SBR, CBR, YBR, V, MU)
%   returns 4 matrices containing the partial derivatives w.r.t. voltage
%   angle and magnitude of the product of a vector MU with the 1st partial
%   derivatives of the square of the magnitude of branch complex power flows.
%   Takes sparse first derivative matrices of complex flow, complex flow
%   vector, sparse connection matrix CBR, sparse branch admittance matrix YBR,
%   voltage vector V and nl x 1 vector of multipliers MU. Output matrices
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
%             d2ASbr_dV2_P(dSbr_dVa, dSbr_dVm, Sbr, Cbr, Ybr, V, mu);
%
%   Here the output matrices correspond to:
%     Haa = (d/dVa (dASbr_dVa.')) * mu
%     Hav = (d/dVm (dASbr_dVa.')) * mu
%     Hva = (d/dVa (dASbr_dVm.')) * mu
%     Hvv = (d/dVm (dASbr_dVm.')) * mu
%
%   See also DSBR_DV_P.
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
nl = length(mu);

diagmu = sparse(1:nl, 1:nl, mu, nl, nl);
diagSbr_conj = sparse(1:nl, 1:nl, conj(Sbr), nl, nl);

[Saa, Sav, Sva, Svv] = d2Sbr_dV2_P(Cbr, Ybr, V, diagSbr_conj * mu);
Haa = 2 * real( Saa + dSbr_dVa.' * diagmu * conj(dSbr_dVa) );
Hva = 2 * real( Sva + dSbr_dVm.' * diagmu * conj(dSbr_dVa) );
Hav = 2 * real( Sav + dSbr_dVa.' * diagmu * conj(dSbr_dVm) );
Hvv = 2 * real( Svv + dSbr_dVm.' * diagmu * conj(dSbr_dVm) );
