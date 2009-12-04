function [Gaa, Gav, Gva, Gvv] = ...
    d2ASbr_dV2(dSbr_dVa, dSbr_dVm, Sbr, Cbr, Ybr, V, lam)
%D2ASBR_DV2   Computes 2nd derivatives of |complex power flow|^2 w.r.t. V.
%   [Gaa, Gav, Gva, Gvv] = d2ASbr_dV2(dSbr_dVa, dSbr_dVm, Sbr, Ybr, V, lam)
%   returns 4 matrices containing the partial derivatives w.r.t. voltage
%   angle and magnitude of the product of a vector lam with the 1st partial
%   derivatives of the square of the magnitude of branch complex power flows.
%   Takes sparse first derivative matrices of complex flow, complex flow
%   vector, sparse connection matrix Cbr, sparse, branch admittance matrix
%   Ybr, voltage vector V and nl x 1 vector lam. Output matrices are sparse.
%
%   e.g f = branch(:, F_BUS);
%     Cf =  sparse(1:nl, f, ones(nl, 1), nl, nb);
%     [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%     [dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = dSbr_dV(branch, Yf, Yt, V)
%     Cbr = Cf;
%     Ybr = Yf;
%     dSbr_dVa = dSf_dVa;
%     dSbr_dVm = dSf_dVm;
%     Sbr = Sf;
%
%   Gaa = (d/dVa (dASbr_dVa.')) * lam
%   Gav = (d/dVm (dASbr_dVa.')) * lam
%   Gva = (d/dVa (dASbr_dVm.')) * lam
%   Gvv = (d/dVm (dASbr_dVm.')) * lam

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define
nb = length(V);
nl = length(lam);

diaglam = sparse(1:nl, 1:nl, lam, nl, nl);
diagSbr_conj = sparse(1:nl, 1:nl, conj(Sbr), nl, nl);

[Saa, Sav, Sva, Svv] = d2Sbr_dV2(Cbr, Ybr, V, diagSbr_conj * lam);
Gaa = 2 * real( Saa + dSbr_dVa.' * diaglam * conj(dSbr_dVa) );
Gva = 2 * real( Sva + dSbr_dVm.' * diaglam * conj(dSbr_dVa) );
Gav = 2 * real( Sav + dSbr_dVa.' * diaglam * conj(dSbr_dVm) );
Gvv = 2 * real( Svv + dSbr_dVm.' * diaglam * conj(dSbr_dVm) );
