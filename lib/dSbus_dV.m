function [dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V)
%DSBUS_DV   Computes partial derivatives of power injection w.r.t. voltage.
%   [dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V) returns two matrices containing
%   partial derivatives of the complex bus power injections w.r.t voltage
%   magnitude and voltage angle respectively (for all buses). If Ybus is a
%   sparse matrix, the return values will be also. The following explains
%	the expressions used to form the matrices:
%
%   S = diag(V) * conj(Ibus) = diag(conj(Ibus)) * V
%
%   Partials of V & Ibus w.r.t. voltage magnitudes
%       dV/dVm = diag(V./abs(V))
%       dI/dVm = Ybus * dV/dVm = Ybus * diag(V./abs(V))
%
%   Partials of V & Ibus w.r.t. voltage angles
%       dV/dVa = j * diag(V)
%       dI/dVa = Ybus * dV/dVa = Ybus * j * diag(V)
%
%   Partials of S w.r.t. voltage magnitudes
%       dS/dVm = diag(V) * conj(dI/dVm) + diag(conj(Ibus)) * dV/dVm
%              = diag(V) * conj(Ybus * diag(V./abs(V)))
%                                       + conj(diag(Ibus)) * diag(V./abs(V))
%
%   Partials of S w.r.t. voltage angles
%       dS/dVa = diag(V) * conj(dI/dVa) + diag(conj(Ibus)) * dV/dVa
%              = diag(V) * conj(Ybus * j * diag(V))
%                                       + conj(diag(Ibus)) * j * diag(V)
%              = -j * diag(V) * conj(Ybus * diag(V))
%                                       + conj(diag(Ibus)) * j * diag(V)
%              = j * diag(V) * conj(diag(Ibus) - Ybus * diag(V))
%

%   MATPOWER Version 2.0
%   by Ray Zimmerman, PSERC Cornell    9/19/97
%   Copyright (c) 1996, 1997 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

j = sqrt(-1);
n = length(V);
Ibus = Ybus * V;

if issparse(Ybus)			%% sparse version (if Ybus is sparse)
	diagV		= spdiags(V, 0, n, n);
	diagIbus	= spdiags(Ibus, 0, n, n);
	diagVnorm	= spdiags(V./abs(V), 0, n, n);
else						%% dense version
	diagV		= diag(V);
	diagIbus	= diag(Ibus);
	diagVnorm	= diag(V./abs(V));
end

dSbus_dVm = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm;
dSbus_dVa = j * diagV * conj(diagIbus - Ybus * diagV);

return;
