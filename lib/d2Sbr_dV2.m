function [Gaa, Gav, Gva, Gvv] = d2Sbr_dV2(Cbr, Ybr, V, lam)
%D2SBR_DV2   Computes 2nd derivatives of complex power flow w.r.t. voltage.
%   [Gaa, Gav, Gva, Gvv] = d2Sbr_dV2(Cbr, Ybr, V, lam) returns 4 matrices
%   containing the partial derivatives w.r.t. voltage angle and magnitude
%   of the product of a vector lam with the 1st partial derivatives of the
%   complex branch power flows. Takes sparse connection matrix Cbr, sparse
%   branch admittance matrix Ybr, voltage vector V and nl x 1 vector lam.
%   Output matrices are sparse.
%
%   e.g f = branch(:, F_BUS);
%       Cf =  sparse(1:nl, f, ones(nl, 1), nl, nb);
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       Cbr = Cf;
%       Ybr = Yf;
%
%   Gaa = (d/dVa (dSbr_dVa.')) * lam
%   Gav = (d/dVm (dSbr_dVa.')) * lam
%   Gva = (d/dVa (dSbr_dVm.')) * lam
%   Gvv = (d/dVm (dSbr_dVm.')) * lam

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define
j = sqrt(-1);
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

Gaa = F - D - E;
Gva = j * G * (B - B.' - D + E);
Gav = Gva.';
Gvv = G * F * G;

return;
