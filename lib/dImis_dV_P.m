function [dImis_dVm, dImis_dVa] = dImis_dV_P(Sbus, Ybus, V)
%dImis_dV_P   Computes partial derivatives of current balance w.r.t. voltage.
%   [DIMIS_DVM, DIMIS_DVA] = dImis_dV_P(SBUS, YBUS, V) returns two matrices containing
%   partial derivatives of the current balance w.r.t voltage
%   magnitude and voltage angle respectively (for all buses). If YBUS is a
%   sparse matrix, the return values will be also. The following explains
%   the expressions used to form the matrices:
%
%   Imis = Ibus - Ispe = Ybus * V - conj(Sbus./V)
%
%   Partials of V & Ibus w.r.t. voltage magnitudes
%       dV/dVm = diag(V./abs(V))
%
%   Partials of V & Ibus w.r.t. voltage angles
%       dV/dVa = j * diag(V)
%
%   Partials of Imis w.r.t. voltage magnitudes
%       dImis/dVm = Ybus * diag(V./abs(V)) + conj( diag( S_sp./(V*abs(V)) ) )
%
%   Partials of Imis w.r.t. voltage angles
%       dImis/dVa = j * (Ybus* diag(V) - conj(diag(S_sp./V)))
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dImis_dVm, dImis_dVa] = dImis_dV_P(Sbus, Ybus, V);
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:

n = length(V);

if issparse(Ybus)
    diagV            = sparse(1:n, 1:n, V, n, n);
    diag_dI_spec_dVm = sparse(1:n, 1:n, Sbus./(V.*abs(V)), n, n);
    diag_dI_spec_dVa = sparse(1:n, 1:n, Sbus./V, n, n); 
    diagVnorm        = sparse(1:n, 1:n, V./abs(V), n, n);
else
    diagV            = diag(V);
    diag_dI_spec_dVm = diag(Sbus./(V.*abs(V)));
    diag_dI_spec_dVa = diag(Sbus./V);
    diagVnorm        = diag(V./abs(V)); 
end

dImis_dVm = Ybus*diagVnorm + conj(diag_dI_spec_dVm);
dImis_dVa = 1j*(Ybus*diagV - conj(diag_dI_spec_dVa));

end