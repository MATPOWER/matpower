function [Haa, Hav, Hva, Hvv] = d2Sbus_dV2(Ybus, V, lam)
%D2SBUS_DV2   Computes 2nd derivatives of power injection w.r.t. voltage.
%   [Haa, Hav, Hva, Hvv] = d2Sbus_dV2(Ybus, V, lam) returns 4 matrices
%   containing the partial derivatives w.r.t. voltage angle and magnitude
%   of the product of a vector lam with the 1st partial derivatives of the
%   complex bus power injections. Takes sparse bus admittance matrix Ybr,
%   voltage vector V and nb x 1 vector lam. Output matrices are sparse.
%
%   Haa = (d/dVa (dSbus_dVa.')) * lam
%   Hav = (d/dVm (dSbus_dVa.')) * lam
%   Hva = (d/dVa (dSbus_dVm.')) * lam
%   Hvv = (d/dVm (dSbus_dVm.')) * lam

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

j = sqrt(-1);
n = length(V);
Ibus    = Ybus * V;
diaglam = sparse(1:n, 1:n, lam, n, n);
diagV   = sparse(1:n, 1:n, V, n, n);

A = sparse(1:n, 1:n, lam .* V, n, n);
B = Ybus * diagV;
C = A * conj(B);
D = Ybus' * diagV;
E = conj(diagV) * (D * diaglam - sparse(1:n, 1:n, D*lam, n, n));
F = C - A * sparse(1:n, 1:n, conj(Ibus), n, n);
G = sparse(1:n, 1:n, ones(n, 1)./abs(V), n, n);

Haa = E + F;
Hva = j * G * (E - F);
Hav = Hva.';
Hvv = G * (C + C.') * G;
