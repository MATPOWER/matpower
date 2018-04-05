function [dSbus_dV1, dSbus_dV2] = dSbus_dV(Ybus, V, vcart)
%DSBUS_DV   Computes partial derivatives of power injection w.r.t. voltage.
%
%   The derivatives can be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 3rd argument.
%
%   [DSBUS_DVA, DSBUS_DVM] = DSBUS_DV(YBUS, V)
%   [DSBUS_DVA, DSBUS_DVM] = DSBUS_DV(YBUS, V, 0)
%
%   Returns two matrices containing partial derivatives of the complex bus
%   power injections w.r.t voltage angle and voltage magnitude, respectively
%   (for all buses).
%
%   [DSBUS_DVR, DSBUS_DVI] = DSBUS_DV(YBUS, V, 1)
%
%   Returns two matrices containing partial derivatives of the complex bus
%   power injections w.r.t the real and imaginary parts of voltage,
%   respectively (for all buses).
%
%   If YBUS is a sparse matrix, the return values will be also. The following
%   explains the expressions used to form the matrices:
%
%   S = diag(V) * conj(Ibus) = diag(conj(Ibus)) * V
%
%   Polar coordinates:
%     Partials of V & Ibus w.r.t. voltage magnitudes
%       dV/dVm = diag(V./abs(V))
%       dI/dVm = Ybus * dV/dVm = Ybus * diag(V./abs(V))
%
%     Partials of V & Ibus w.r.t. voltage angles
%       dV/dVa = j * diag(V)
%       dI/dVa = Ybus * dV/dVa = Ybus * j * diag(V)
%
%     Partials of S w.r.t. voltage magnitudes
%       dS/dVm = diag(V) * conj(dI/dVm) + diag(conj(Ibus)) * dV/dVm
%              = diag(V) * conj(Ybus * diag(V./abs(V)))
%                                       + conj(diag(Ibus)) * diag(V./abs(V))
%
%     Partials of S w.r.t. voltage angles
%       dS/dVa = diag(V) * conj(dI/dVa) + diag(conj(Ibus)) * dV/dVa
%              = diag(V) * conj(Ybus * j * diag(V))
%                                       + conj(diag(Ibus)) * j * diag(V)
%              = -j * diag(V) * conj(Ybus * diag(V))
%                                       + conj(diag(Ibus)) * j * diag(V)
%              = j * diag(V) * conj(diag(Ibus) - Ybus * diag(V))
%
%   Cartesian coordinates:
%     Partials of V & Ibus w.r.t. real part of complex voltage
%       dV/dVr = diag(ones(n,1))
%       dI/dVr = Ybus * dV/dVr = Ybus
%
%     Partials of V & Ibus w.r.t. imaginary part of complex voltage
%       dV/dVi = j * diag(ones(n,1))
%       dI/dVi = Ybus * dV/dVi = Ybus * j
%
%     Partials of S w.r.t. real part of complex voltage
%       dS/dVr = diag(V) * conj(dI/dVr) + diag(conj(Ibus)) * dV/dVr
%              = diag(V) * conj(Ybus) + conj(diag(Ibus))
%
%     Partials of S w.r.t. imaginary part of complex voltage
%       dS/dVi = diag(V) * conj(dI/dVi) + diag(conj(Ibus)) * dV/dVi
%              = j * (conj(diag(Ibus)) - diag(V) conj(Ybus))
%
%   Examples:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dSbus_dVa, dSbus_dVm] = dSbus_dV(Ybus, V);
%       [dSbus_dVr, dSbus_dVi] = dSbus_dV(Ybus, V, 1);
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010.
%             http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf
%   [TN4]  B. Sereeter and R. D. Zimmerman, "AC Power Flows and their
%          Derivatives using Complex Matrix Notation and Cartesian
%          Coordinate Voltages," MATPOWER Technical Note 4, April 2018.
%             http://www.pserc.cornell.edu/matpower/
%                                           TN4-OPF-Derivatives-Cartesian.pdf

%   MATPOWER
%   Copyright (c) 1996-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default input args
if nargin < 3
    vcart = 0;      %% default to polar coordinates
end

n = length(V);
Ibus = Ybus * V;

if issparse(Ybus)           %% sparse version (if Ybus is sparse)
    diagV       = sparse(1:n, 1:n, V, n, n);
    diagIbus    = sparse(1:n, 1:n, Ibus, n, n);
    if ~vcart
        diagVnorm   = sparse(1:n, 1:n, V./abs(V), n, n);
    end
else                        %% dense version
    diagV       = diag(V);
    diagIbus    = diag(Ibus);
    if ~vcart
        diagVnorm   = diag(V./abs(V));
    end
end

if vcart
    dSbus_dV1 = conj(diagIbus) + diagV * conj(Ybus);        %% dSbus/dVr
    dSbus_dV2 = 1j * (conj(diagIbus) - diagV * conj(Ybus)); %% dSbus/dVi
else
    dSbus_dV1 = 1j * diagV * conj(diagIbus - Ybus * diagV);                     %% dSbus/dVa
    dSbus_dV2 = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm;    %% dSbus/dVm
end
