function [Gii, Gir, Gri, Grr] = d2Imis_dV2_C(Sbus, Ybus, V, lam)
%d2Imis_dV2_C   Computes 2nd derivatives of current balance w.r.t. voltage.
%   [GII, GIR, GRI, GRR] = d2Imis_dV2_C(YBUS, V, LAM) returns 4 matrices
%   containing the partial derivatives w.r.t. real and imaginary part of complex
%   of the product of a vector LAM with the 1st partial derivatives of the
%   complex current balance. Takes sparse bus admittance matrix YBUS,
%   voltage vector V and nb x 1 vector of multipliers LAM. Output matrices
%   are sparse.
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [Gii, Gir, Gri, Grr] = d2Imis_dV2_C(Ybus, V, lam);
%
%   Here the output matrices correspond to:
%       Gii = (d/dVi (dImis_dVi.')) * lam
%       Gir = (d/dVr (dImis_dVi.')) * lam
%       Gri = (d/dVi (dImis_dVr.')) * lam
%       Grr = (d/dVr (dImis_dVr.')) * lam
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER, see:

nb = length(V);
diagLamS  = sparse(1:nb, 1:nb, Sbus.*lam, nb, nb);    
diagV3    = sparse(1:nb, 1:nb, 1./(V.^3), nb, nb);    
    
Gii   = 2*conj(diagLamS*diagV3);
Grr   = - G_Vi2;
Gir = 1j*G_Vi2;
Gri = Gir;
end
