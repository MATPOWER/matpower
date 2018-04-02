function [H11, H12, H21, H22] = d2Ibr_dV2(Ybr, V, mu, vcart)
%D2IBR_DV2   Computes 2nd derivatives of complex branch current w.r.t. voltage.
%
%   The derivatives can be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 4th argument.
%
%   [HAA, HAV, HVA, HVV] = D2IBR_DV2(YBR, V, MU)
%   [HAA, HAV, HVA, HVV] = D2IBR_DV2(YBR, V, MU, 0)
%
%   Returns 4 matrices containing the partial derivatives w.r.t. voltage angle
%   and magnitude of the product of a vector MU with the 1st partial
%   derivatives of the complex branch currents.
%
%   [HRR, HRI, HIR, HII] = D2IBR_DV2(YBR, V, MU, 1)
%
%   Returns 4 matrices (all zeros) containing the partial derivatives w.r.t.
%   real and imaginary part of complex voltage of the product of a vector MU
%   with the 1st partial derivatives of the complex branch currents.
%
%   Takes sparse branch admittance matrix YBR, voltage vector V and nl x 1
%   vector of multipliers MU. Output matrices are sparse.
%
%   Examples:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       Ybr = Yf;
%       [Haa, Hav, Hva, Hvv] = d2Ibr_dV2(Ybr, V, mu);
%
%   Here the output matrices correspond to:
%       Haa = d/dVa (dIbr_dVa.' * mu)
%       Hav = d/dVm (dIbr_dVa.' * mu)
%       Hva = d/dVa (dIbr_dVm.' * mu)
%       Hvv = d/dVm (dIbr_dVm.' * mu)
%
%       [Hrr, Hri, Hir, Hii] = d2Ibr_dV2(Ybr, V, mu, 1);
%
%   Here the output matrices correspond to:
%       Hrr = d/dVr (dIbr_dVr.' * mu)
%       Hri = d/dVi (dIbr_dVr.' * mu)
%       Hir = d/dVr (dIbr_dVi.' * mu)
%       Hii = d/dVi (dIbr_dVi.' * mu)
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
%   Copyright (c) 2008-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default input args
if nargin < 4
    vcart = 0;      %% default to polar coordinates
end

%% define
nb = length(V);

if vcart
    H11 = sparse(nb, nb);
    H12 = H11;
    H21 = H11;
    H22 = H11;
else
    diagInvVm = sparse(1:nb, 1:nb, ones(nb, 1)./abs(V), nb, nb);

    H11 = sparse(1:nb, 1:nb, -(Ybr.' * mu) .* V, nb, nb);
    H21 = -1j * H11 * diagInvVm;
    H12 = H21;
    H22 = sparse(nb, nb);
end
