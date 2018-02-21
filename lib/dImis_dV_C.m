function [dImis_dVr, dImis_dVi] = dImis_dV_C(Sbus, Ybus, V)
%dImis_dV_C   Computes partial derivatives of current balance w.r.t. voltage.
%   [DSMIS_DVR, DIMIS_DVI] = dImis_dV_C(SBUS, YBUS, V) returns two matrices containing
%   partial derivatives of the current balance w.r.t real and 
%   imaginary part of complex voltage respectively (for all buses). If YBUS is a
%   sparse matrix, the return values will be also. The following explains
%   the expressions used to form the matrices:
%
%   Imis = Ibus - Ispe = Ybus * V - conj(Sbus./V)
%
%   Partials of V & Ibus w.r.t. real part of complex voltage
%       dV/dVr = diag(ones(n,1))
%
%   Partials of V & Ibus w.r.t. imaginary part of complex voltage 
%       dV/dVi = j * diag(ones(n,1))
%
%   Partials of S w.r.t. real part of complex voltage
%       dImis/dVr = Ybus + conj(diag(Sbus./(V.^2)))
%
%   Partials of S w.r.t. imaginary part of complex voltage 
%       dImis/dVi = j * (Ybus - conj(diag(Sbus./(V.^2))) )
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dImis_dVr, dImis_dVi] = dImis_dV_C(Sbus, Ybus, V);
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:

n = length(V);

if issparse(Ybus)
   diagSV2     = sparse(1:n, 1:n, Sbus./(V.^2), n, n);
else
   diagSV2     = diag(Sbus./(V.^2));
end

dImis_dVr = Ybus + conj(diagSV2);
dImis_dVi = 1j*(Ybus - conj(diagSV2)); 

end