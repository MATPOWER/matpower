function [G11, G12, G21, G22] = d2Sbus_dV2(Ybus, V, lam, vcart)
%D2SBUS_DV2   Computes 2nd derivatives of power injection w.r.t. voltage.
%
%   The derivatives can be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 4th argument.
%
%   [GAA, GAV, GVA, GVV] = D2SBUS_DV2(YBUS, V, LAM)
%   [GAA, GAV, GVA, GVV] = D2SBUS_DV2(YBUS, V, LAM, 0)
%
%   Returns 4 matrices containing the partial derivatives w.r.t. voltage angle
%   and magnitude of the product of a vector LAM with the 1st partial
%   derivatives of the complex bus power injections.
%
%   [GRR, GIR, GIR, GII] = D2SBUS_DV2(YBUS, V, LAM, 1)
%
%   Returns 4 matrices containing the partial derivatives w.r.t. real and
%   imaginary parts of voltage of the product of a vector LAM with the 1st
%   partial derivatives of the complex bus power injections.
%
%   Takes sparse bus admittance matrix YBUS, voltage vector V and nb x 1 vector
%   of multipliers LAM. Output matrices are sparse.
%
%   Examples:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [Gaa, Gav, Gva, Gvv] = d2Sbus_dV2(Ybus, V, lam);
%
%       Here the output matrices correspond to:
%           Gaa = d/dVa (dSbus_dVa.' * lam)
%           Gav = d/dVm (dSbus_dVa.' * lam)
%           Gva = d/dVa (dSbus_dVm.' * lam)
%           Gvv = d/dVm (dSbus_dVm.' * lam)
%
%       [Grr, Gri, Gir, Gii] = d2Sbus_dV2(Ybus, V, lam, 1);
%
%       Here the output matrices correspond to:
%           Grr = d/dVr (dSbus_dVr.' * lam)
%           Gri = d/dVi (dSbus_dVr.' * lam)
%           Gir = d/dVr (dSbus_dVi.' * lam)
%           Gii = d/dVi (dSbus_dVi.' * lam)
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER, see:
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
%   Copyright (c) 2008-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default input args
if nargin < 4
    vcart = 0;      %% default to polar coordinates
end

n = length(V);

diaglam = sparse(1:n, 1:n, lam, n, n);
if vcart
    E = diaglam * conj(Ybus);
    F = E + E.';
    G = 1j * (E - E.');

    G11 = F;        %% Grr
    G21 = G;        %% Gir
    G12 = G21.';    %% Gri
    G22 = G11;      %% Gii
else
    Ibus    = Ybus * V;
    diagV   = sparse(1:n, 1:n, V, n, n);

    A = sparse(1:n, 1:n, lam .* V, n, n);
    B = Ybus * diagV;
    C = A * conj(B);
    D = Ybus' * diagV;
    E = conj(diagV) * (D * diaglam - sparse(1:n, 1:n, D*lam, n, n));
    F = C - A * sparse(1:n, 1:n, conj(Ibus), n, n);
    G = sparse(1:n, 1:n, ones(n, 1)./abs(V), n, n);

    G11 = E + F;
    G21 = 1j * G * (E - F);
    G12 = G21.';
    G22 = G * (C + C.') * G;
end
