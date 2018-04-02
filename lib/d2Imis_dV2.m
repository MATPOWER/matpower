function [G11, G12, G21, G22] = d2Imis_dV2(Sbus, Ybus, V, lam, vcart)
%D2IMIS_DV2   Computes 2nd derivatives of current balance w.r.t. voltage.
%
%   The derivatives can be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 5th argument.
%
%   [GAA, GAV, GVA, GVV] = D2IMIS_DV2(SBUS, YBUS, V, LAM)
%   [GAA, GAV, GVA, GVV] = D2IMIS_DV2(SBUS, YBUS, V, LAM, 0)
%
%   Returns 4 matrices containing the partial derivatives w.r.t. voltage angle
%   and magnitude of the product of a vector LAM with the 1st partial
%   derivatives of the complex bus current balance.
%
%   [GRR, GIR, GIR, GII] = D2IMIS_DV2(SBUS, YBUS, V, LAM, 1)
%
%   Returns 4 matrices containing the partial derivatives w.r.t. real and
%   imaginary parts of voltage of the product of a vector LAM with the 1st
%   partial derivatives of the complex bus current balance.
%
%   Takes bus complex power injection (gen-load) vector, sparse bus admittance
%   matrix YBUS, voltage vector V and nb x 1 vector of multipliers LAM. Output
%   matrices are sparse.
%
%   Examples:
%       Sbus = makeSbus(baseMVA, bus, gen);
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [Gaa, Gav, Gva, Gvv] = d2Imis_dV2(Sbus, Ybus, V, lam);
%
%       Here the output matrices correspond to:
%           Gaa = d/dVa (dImis_dVa.' * lam)
%           Gav = d/dVm (dImis_dVa.' * lam)
%           Gva = d/dVa (dImis_dVm.' * lam)
%           Gvv = d/dVm (dImis_dVm.' * lam)
%
%       [Grr, Gri, Gir, Gii] = d2Imis_dV2(Sbus, Ybus, V, lam, 1);
%
%       Here the output matrices correspond to:
%           Grr = d/dVr (dImis_dVr.' * lam)
%           Gri = d/dVi (dImis_dVr.' * lam)
%           Gir = d/dVr (dImis_dVi.' * lam)
%           Gii = d/dVi (dImis_dVi.' * lam)
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER, see:
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
if nargin < 5
    vcart = 0;      %% default to polar coordinates
end

nb = length(V);

if vcart
    C = 2 * lam .* conj(Sbus./(V.^3));

    G22 = sparse(1:nb, 1:nb, C, nb, nb);
    G11 = -G22;
    G12 = 1j*G22;
    G21 = G12;
else
    absV = abs(V);
    diagV    = sparse(1:nb, 1:nb, V, nb, nb);
    diagV1   = sparse(1:nb, 1:nb, 1./V, nb, nb);
    diagVmV1 = sparse(1:nb, 1:nb, 1./(V.*absV), nb, nb);
    diagVmV2 = sparse(1:nb, 1:nb, 1./(V.*(absV.^2)), nb, nb);
    diagLamS = sparse(1:nb, 1:nb, lam.*Sbus, nb, nb);
    diagYlam = sparse(1:nb, 1:nb, Ybus.'*lam, nb, nb);
    diagE    = sparse(1:nb, 1:nb, V./absV, nb, nb);

    G11   = - diagYlam* diagV + conj(diagLamS*diagV1);
    G22   = - 2 * conj(diagLamS*diagVmV2);
    G21 = 1j * ( diagYlam* diagE + conj(diagLamS*diagVmV1)) ;
    G12 = G21.';
end
