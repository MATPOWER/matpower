function [dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V)
%DSBUS_DV   Computes partial derivatives of power injection w.r.t. voltage.
%   [DSBUS_DVM, DSBUS_DVA] = DSBUS_DV(YBUS, V) returns two matrices containing
%   partial derivatives of the complex bus power injections w.r.t voltage
%   magnitude and voltage angle respectively (for all buses). If YBUS is a
%   sparse matrix, the return values will be also. The following explains
%   the expressions used to form the matrices:
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
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V);
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010.
%             http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

n = length(V);
Ibus = Ybus * V;

if issparse(Ybus)           %% sparse version (if Ybus is sparse)
    diagV       = sparse(1:n, 1:n, V, n, n);
    diagIbus    = sparse(1:n, 1:n, Ibus, n, n);
    diagVnorm   = sparse(1:n, 1:n, V./abs(V), n, n);
else                        %% dense version
    diagV       = diag(V);
    diagIbus    = diag(Ibus);
    diagVnorm   = diag(V./abs(V));
end

dSbus_dVm = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm;
dSbus_dVa = 1j * diagV * conj(diagIbus - Ybus * diagV);
