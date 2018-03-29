function [Grr, Gri, Gir, Gii] = d2Imis_dV2_C(Sbus, Ybus, V, lam)
%d2Imis_dV2_C   Computes 2nd derivatives of current balance w.r.t. voltage.
%   [GRR, GRI, GIR, GII] = d2Imis_dV2_C(YBUS, V, LAM) returns 4 matrices
%   containing the partial derivatives w.r.t. real and imaginary part of complex
%   of the product of a vector LAM with the 1st partial derivatives of the
%   complex current balance. Takes sparse bus admittance matrix YBUS,
%   voltage vector V and nb x 1 vector of multipliers LAM. Output matrices
%   are sparse.
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [Grr, Gri, Gir, Gii] = d2Imis_dV2_C(Ybus, V, lam);
%
%   Here the output matrices correspond to:
%       Grr = (d/dVr (dImis_dVr.')) * lam
%       Gri = (d/dVi (dImis_dVr.')) * lam
%       Gir = (d/dVr (dImis_dVi.')) * lam
%       Gii = (d/dVi (dImis_dVi.')) * lam
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER, see:

nb = length(V);

C = 2 * lam .* conj(Sbus./(V.^3));

Gii = sparse(1:nb, 1:nb, C, nb, nb);
Grr = -Gii;
Gri = 1j*Gii;
Gir = Gri;
