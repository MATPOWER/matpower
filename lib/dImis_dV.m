function [dImis_dV1, dImis_dV2] = dImis_dV(Sbus, Ybus, V, vcart)
%DIMIS_DV   Computes partial derivatives of current balance w.r.t. voltage.
%
%   The derivatives can be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 3rd argument.
%
%   [DIMIS_DVM, DIMIS_DVA] = DIMIS_DV(SBUS, YBUS, V)
%   [DIMIS_DVM, DIMIS_DVA] = DIMIS_DV(SBUS, YBUS, V, 0)
%
%   Returns two matrices containing partial derivatives of the complex bus
%   current balance w.r.t voltage magnitude and voltage angle, respectively
%   (for all buses).
%
%   [DIMIS_DVR, DIMIS_DVI] = DIMIS_DV(SBUS, YBUS, V, 1)
%
%   Returns two matrices containing partial derivatives of the complex bus
%   current balance w.r.t the real and imaginary parts of voltage,
%   respectively (for all buses).
%
%   If YBUS is a sparse matrix, the return values will be also. The following
%   explains the expressions used to form the matrices:
%
%   Imis = Ibus + Idg = Ybus * V - conj(Sbus./V)
%
%   Polar coordinates:
%     Partials of V & Ibus w.r.t. voltage angles
%       dV/dVa = j * diag(V)
%       dI/dVa = Ybus * dV/dVa = Ybus * j * diag(V)
%
%     Partials of V & Ibus w.r.t. voltage magnitudes
%       dV/dVm = diag(V./abs(V))
%       dI/dVm = Ybus * dV/dVm = Ybus * diag(V./abs(V))
%
%     Partials of Imis w.r.t. voltage angles
%       dImis/dVa = j * (Ybus* diag(V) - diag(conj(Sbus./V)))
%
%     Partials of Imis w.r.t. voltage magnitudes
%       dImis/dVm = Ybus * diag(V./abs(V)) +  diag(conj(Sbus./(V*abs(V))))
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
%     Partials of Imis w.r.t. real part of complex voltage
%       dImis/dVr = Ybus + conj(diag(Sbus./(V.^2)))
%
%     Partials of S w.r.t. imaginary part of complex voltage
%       dImis/dVi = j * (Ybus - diag(conj(Sbus./(V.^2))))
%
%   Examples:
%       Sbus = makeSbus(baseMVA, bus, gen);
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dImis_dVm, dImis_dVa] = dImis_dV(Sbus, Ybus, V);
%       [dImis_dVr, dImis_dVi] = dImis_dV(Sbus, Ybus, V, 1);
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010.
%             http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf
%   [TN3]  B. Sereeter and R. D. Zimmerman, "Addendum to AC Power Flows and
%          their Derivatives using Complex Matrix Notation: Nodal Current
%          Balance," MATPOWER Technical Note 3, April 2018.
%             http://www.pserc. cornell.edu/matpower/
%                                           TN3-More-OPF-Derivatives.pdf
%   [TN4]  B. Sereeter and R. D. Zimmerman, "AC Power Flows and their
%          Derivatives using Complex Matrix Notation and Cartesian
%          Coordinate Voltages," MATPOWER Technical Note 4, April 2018.
%             http://www.pserc.cornell.edu/matpower/
%                                           TN4-OPF-Derivatives-Cartesian.pdf

%   MATPOWER
%   Copyright (c) 2018, Power Systems Engineering Research Center (PSERC)
%   by Baljinnyam Sereeter, Delft University of Technology
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default input args
if nargin < 4
    vcart = 0;      %% default to polar coordinates
end

n = length(V);

if vcart
    if issparse(Ybus)
       diagSV2c = sparse(1:n, 1:n, conj(Sbus./(V.^2)), n, n);
    else
       diagSV2c = diag(conj(Sbus./(V.^2)));
    end

    dImis_dV1 = Ybus + diagSV2c;            %% dImis/dVr
    dImis_dV2 = 1j*(Ybus - diagSV2c);       %% dImis/dVi
else
    Vm = abs(V);
    Vnorm = V./Vm;
    Ibus = conj(Sbus./V);
    if issparse(Ybus)
        diagV       = sparse(1:n, 1:n, V, n, n);
        diagIbus    = sparse(1:n, 1:n, Ibus, n, n);
        diagIbusVm  = sparse(1:n, 1:n, Ibus./Vm, n, n);
        diagVnorm   = sparse(1:n, 1:n, V./abs(V), n, n);
    else
        diagV       = diag(V);
        diagIbus    = diag(Ibus);
        diagIbusVm  = diag(Ibus./Vm);
        diagVnorm   = diag(V./abs(V));
    end

    dImis_dV1 = 1j*(Ybus*diagV - diagIbus);     %% dImis/dVa
    dImis_dV2 = Ybus*diagVnorm + diagIbusVm;    %% dImis/dVm
end
