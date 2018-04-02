function [H11, H12, H21, H22] = d2Sbr_dV2(Cbr, Ybr, V, mu, vcart)
%D2SBR_DV2   Computes 2nd derivatives of complex brch power flow w.r.t. voltage.
%
%   The derivatives can be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 5th argument.
%
%   [HAA, HAV, HVA, HVV] = D2SBR_DV2(CBR, YBR, V, MU)
%   [HAA, HAV, HVA, HVV] = D2SBR_DV2(CBR, YBR, V, MU, 0)
%
%   Returns 4 matrices containing the partial derivatives w.r.t. voltage angle
%   and magnitude of the product of a vector MU with the 1st partial
%   derivatives of the complex branch power flows.
%
%   [HRR, HRI, HIR, HII] = d2Sbr_dV2(CBR, YBR, V, MU, 1)
%
%   Returns 4 matrices containing the partial derivatives w.r.t. real and
%   imaginary part of complex voltage of the product of a vector MU with the
%   1st partial derivatives of the complex branch power flows.
%
%   Takes sparse connection matrix CBR, sparse branch admittance matrix YBR,
%   voltage vector V and nl x 1 vector of multipliers MU. Output matrices are
%   sparse.
%
%   Examples:
%       f = branch(:, F_BUS);
%       Cf =  sparse(1:nl, f, ones(nl, 1), nl, nb);
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       Cbr = Cf;
%       Ybr = Yf;
%       [Haa, Hav, Hva, Hvv] = d2Sbr_dV2(Cbr, Ybr, V, mu);
%
%       Here the output matrices correspond to:
%           Haa = d/dVa (dSbr_dVa.' * mu)
%           Hav = d/dVm (dSbr_dVa.' * mu)
%           Hva = d/dVa (dSbr_dVm.' * mu)
%           Hvv = d/dVm (dSbr_dVm.' * mu)
%
%       [Hrr, Hri, Hir, Hii] = d2Sbr_dV2(Cbr, Ybr, V, mu, 1);
%
%       Here the output matrices correspond to:
%           Hrr = d/dVr (dSbr_dVr.' * mu)
%           Hri = d/dVi (dSbr_dVr.' * mu)
%           Hir = d/dVr (dSbr_dVi.' * mu)
%           Hii = d/dVi (dSbr_dVi.' * mu)
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
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default input args
if nargin < 5
    vcart = 0;      %% default to polar coordinates
end

nl = length(mu);
nb = length(V);

A = Ybr' * sparse(1:nl, 1:nl, mu, nl, nl) * Cbr;
if vcart
    H11 = A + A.';
    H12 = 1j * (A - A.');
    H21 = -H12;
    H22 = H11;
else
    diagV  = sparse(1:nb, 1:nb, V, nb, nb);

    B = conj(diagV) * A * diagV;
    D = sparse(1:nb, 1:nb, (A*V) .* conj(V), nb, nb);
    E = sparse(1:nb, 1:nb, (A.'*conj(V)) .* V, nb, nb);
    F = B + B.';
    G = sparse(1:nb, 1:nb, ones(nb, 1)./abs(V), nb, nb);

    H11 = F - D - E;
    H21 = 1j * G * (B - B.' - D + E);
    H12 = H21.';
    H22 = G * F * G;
end
