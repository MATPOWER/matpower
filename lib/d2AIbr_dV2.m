function [Gaa, Gav, Gva, Gvv] = ...
    d2AIbr_dV2(dIbr_dVa, dIbr_dVm, Ibr, Ybr, V, lam)
%D2AIBR_DV2   Computes 2nd derivatives of |complex current|^2 w.r.t. V.
%   [Gaa, Gav, Gva, Gvv] = d2AIbr_dV2(dIbr_dVa, dIbr_dVm, Ibr, Ybr, V, lam)
%   returns 4 matrices containing the partial derivatives w.r.t. voltage
%   angle and magnitude of the product of a vector lam with the 1st partial
%   derivatives of the square of the magnitude of the branch currents.
%   Takes sparse first derivative matrices of complex flow, complex flow
%   vector, sparse connection matrix Cbr, sparse, branch admittance matrix
%   Ybr, voltage vector V and nl x 1 vector lam. Output matrices are sparse.
%
%   e.g f = branch(:, F_BUS);
%     Cf =  sparse(1:nl, f, ones(nl, 1), nl, nb);
%     [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%     [dIf_dVa, dIf_dVm, dIt_dVa, dIt_dVm, If, It] = dIbr_dV(branch, Yf, Yt, V)
%     Cbr = Cf;
%     Ybr = Yf;
%     dIbr_dVa = dIf_dVa;
%     dIbr_dVm = dIf_dVm;
%     Ibr = If;
%
%   Gaa = (d/dVa (dAIbr_dVa.')) * lam
%   Gav = (d/dVm (dAIbr_dVa.')) * lam
%   Gva = (d/dVa (dAIbr_dVm.')) * lam
%   Gvv = (d/dVm (dAIbr_dVm.')) * lam

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define
j = sqrt(-1);
nb = length(V);
nl = length(lam);

diaglam = sparse(1:nl, 1:nl, lam, nl, nl);
diagIbr_conj = sparse(1:nl, 1:nl, conj(Ibr), nl, nl);

[Iaa, Iav, Iva, Ivv] = d2Ibr_dV2(Ybr, V, diagIbr_conj * lam);
Gaa = 2 * real( Iaa + dIbr_dVa.' * diaglam * conj(dIbr_dVa) );
Gva = 2 * real( Iva + dIbr_dVm.' * diaglam * conj(dIbr_dVa) );
Gav = 2 * real( Iav + dIbr_dVa.' * diaglam * conj(dIbr_dVm) );
Gvv = 2 * real( Ivv + dIbr_dVm.' * diaglam * conj(dIbr_dVm) );

return;
