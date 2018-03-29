function [Grr, Gri, Gir, Gii] = d2Sbus_dV2_C(Ybus, V, lam)
%d2Sbus_dV2_C   Computes 2nd derivatives of power injection w.r.t. voltage.
%   [GRR, GIR, GIR, GII] = d2Sbus_dV2_C(YBUS, V, LAM) returns 4 matrices
%   containing the partial derivatives w.r.t. real and imaginary part of complex
%   of the product of a vector LAM with the 1st partial derivatives of the
%   complex bus power injections. Takes sparse bus admittance matrix YBUS,
%   voltage vector V and nb x 1 vector of multipliers LAM. Output matrices
%   are sparse.
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [Grr, Gri, Gir, Gii] = d2Sbus_dV2_C(Ybus, V, lam);
%
%   Here the output matrices correspond to:
%       Grr = (d/dVr (dSbus_dVr.')) * lam
%       Gri = (d/dVi (dSbus_dVr.')) * lam
%       Gir = (d/dVr (dSbus_dVi.')) * lam
%       Gii = (d/dVi (dSbus_dVi.')) * lam
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER, see:

n = length(V);

diaglam  = sparse(1:n, 1:n, lam, n, n);
LamYbus  = diaglam*conj(Ybus);
YbusTLam = (Ybus')*diaglam;

Grr = LamYbus + YbusTLam;
Gri = 1j*(YbusTLam - LamYbus);
Gir = -Gri;
Gii = Grr;
