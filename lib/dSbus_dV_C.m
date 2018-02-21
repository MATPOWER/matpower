function [dSbus_dVr, dSbus_dVi] = dSbus_dV_C(Ybus, V)
%dSbus_dV_C   Computes partial derivatives of power injection w.r.t. voltage.
%   [DSBUS_DVR, DSBUS_DVI] = dSbus_dV_C(YBUS, V) returns two matrices containing
%   partial derivatives of the complex bus power injections w.r.t real and 
%   imaginary part of complex voltage respectively (for all buses). If YBUS is a
%   sparse matrix, the return values will be also. The following explains
%   the expressions used to form the matrices:
%
%   S = diag(V) * conj(Ibus) = diag(conj(Ibus)) * V
%
%   Partials of V & Ibus w.r.t. real part of complex voltage
%       dV/dVr = diag(ones(n,1))
%       dI/dVr = Ybus * dV/dVr = Ybus 
%
%   Partials of V & Ibus w.r.t. imaginary part of complex voltage 
%       dV/dVi = j * diag(ones(n,1))
%       dI/dVi = Ybus * dV/dVi = Ybus * j 
%
%   Partials of S w.r.t. real part of complex voltage
%       dS/dVr = diag(V) * conj(dI/dVr) + diag(conj(Ibus)) * dV/dVr
%              = diag(V) * conj(Ybus) + conj(diag(Ibus)) 
%
%   Partials of S w.r.t. imaginary part of complex voltage 
%       dS/dVi = diag(V) * conj(dI/dVi) + diag(conj(Ibus)) * dV/dVi
%              = j * (conj(diag(Ibus)) - diag(V) conj(Ybus))
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dSbus_dVr, dSbus_dVi] = dSbus_dV_C(Ybus, V);
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:

n = length(V);
Ibus = Ybus * V;
if issparse(Ybus)           %% sparse version (if Ybus is sparse)
    diagV       = sparse(1:n, 1:n, V, n, n);
    diagIbus    = sparse(1:n, 1:n, Ibus, n, n);
else                        %% dense version
    diagV       = diag(V);
    diagIbus    = diag(Ibus);
end
dSbus_dVi = 1j* (conj(diagIbus) - diagV*conj(Ybus));   
dSbus_dVr = conj(diagIbus) + diagV* conj(Ybus);   
end
