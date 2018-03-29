function [Gaa, Gav, Gva, Gvv] = d2Imis_dV2_P(Sbus, Ybus, V, lam)
%d2Imis_dV2_P   Computes 2nd derivatives of current balance w.r.t. voltage.
%   [GAA, GAV, GVA, GVV] = d2Imis_dV2_P(YBUS, V, LAM) returns 4 matrices
%   containing the partial derivatives w.r.t. voltage angle and magnitude
%   of the product of a vector LAM with the 1st partial derivatives of the
%   complex current balance. Takes sparse bus admittance matrix YBUS,
%   voltage vector V and nb x 1 vector of multipliers LAM. Output matrices
%   are sparse.
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [Gaa, Gav, Gva, Gvv] = d2Imis_dV2_P(Ybus, V, lam);
%
%   Here the output matrices correspond to:
%       Gaa = (d/dVa (dImis_dVa.')) * lam
%       Gav = (d/dVm (dImis_dVa.')) * lam
%       Gva = (d/dVa (dImis_dVm.')) * lam
%       Gvv = (d/dVm (dImis_dVm.')) * lam
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER, see:

nb = length(V);
absV = abs(V);
diagV    = sparse(1:nb, 1:nb, V, nb, nb);
diagV1   = sparse(1:nb, 1:nb, 1./V, nb, nb);
diagVmV1 = sparse(1:nb, 1:nb, 1./(V.*absV), nb, nb);
diagVmV2 = sparse(1:nb, 1:nb, 1./(V.*(absV.^2)), nb, nb);
diagLamS = sparse(1:nb, 1:nb, lam.*Sbus, nb, nb);
diagYlam = sparse(1:nb, 1:nb, Ybus.'*lam, nb, nb);
diagE    = sparse(1:nb, 1:nb, V./absV, nb, nb);

Gaa   = - diagYlam* diagV + conj(diagLamS*diagV1);
Gvv   = - 2 * conj(diagLamS*diagVmV2);
Gva = 1j * ( diagYlam* diagE + conj(diagLamS*diagVmV1)) ;
Gav = Gva.';
