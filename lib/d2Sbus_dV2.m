function [Gaa, Gav, Gva, Gvv] = d2Sbus_dV2(Ybus, V, lam)
%D2SBUS_DV2   Computes 2nd derivatives of power injection w.r.t. voltage.
%   [GAA, GAV, GVA, GVV] = D2SBUS_DV2(YBUS, V, LAM) returns 4 matrices
%   containing the partial derivatives w.r.t. voltage angle and magnitude
%   of the product of a vector LAM with the 1st partial derivatives of the
%   complex bus power injections. Takes sparse bus admittance matrix YBUS,
%   voltage vector V and nb x 1 vector of multipliers LAM. Output matrices
%   are sparse.
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [Gaa, Gav, Gva, Gvv] = d2Sbus_dV2(Ybus, V, lam);
%
%   Here the output matrices correspond to:
%       Gaa = (d/dVa (dSbus_dVa.')) * lam
%       Gav = (d/dVm (dSbus_dVa.')) * lam
%       Gva = (d/dVa (dSbus_dVm.')) * lam
%       Gvv = (d/dVm (dSbus_dVm.')) * lam
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER, see:
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

Gaa = E + F;
Gva = 1j * G * (E - F);
Gav = Gva.';
Gvv = G * (C + C.') * G;
