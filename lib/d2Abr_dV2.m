function [H11, H12, H21, H22] = d2Abr_dV2(d2F_dV2, dF_dV1, dF_dV2, F, V, mu)
%D2ABR_DV2   Computes 2nd derivatives of |branch flow|^2 w.r.t. V.
%
%   The derivatives can be take with respect to polar or cartesian coordinates
%   of voltage, depending on the first 3 arguments. Flows could be complex
%   current or complex or real power. Notation below is based on complex power.
%
%   [H11, H12, H21, H22] = D2ABR_DV2(D2F_DV2, DF_DV1, DF_DV2, F, V, MU)
%
%   Returns 4 matrices containing the partial derivatives w.r.t. voltage
%   components (angle, magnitude or real, imaginary) of the product of a
%   vector MU with the 1st partial derivatives of the square of the magnitude
%   of branch flows.
%
%   Takes as inputs a handle to a function that evaluates the 2nd derivatives
%   of the flows (with args V and mu only), sparse first derivative matrices
%   of flow, flow vector, voltage vector V and nl x 1 vector of multipliers
%   MU. Output matrices are sparse.
%
%   Example:
%       f = branch(:, F_BUS);
%       Cf =  sparse(1:nl, f, ones(nl, 1), nl, nb);
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = ...
%               dSbr_dV(branch, Yf, Yt, V);
%       dF_dV1 = dSf_dV1;
%       dF_dV2 = dSf_dV2;
%       F = Sf;
%       d2F_dV2 = @(V, mu)d2Sbr_dV2(Cf, Yf, V, mu, 0);
%       [H11, H12, H21, H22] = ...
%             d2Abr_dV2(d2F_dV2, dF_dV1, dF_dV2, F, V, mu);
%
%   Here the output matrices correspond to:
%     H11 = d/dV1 (dAF_dV1.' * mu)
%     H12 = d/dV2 (dAF_dV1.' * mu)
%     H21 = d/dV1 (dAF_dV2.' * mu)
%     H22 = d/dV2 (dAF_dV2.' * mu)
%
%   See also DABR_DV, DIBR_DV, DSBR_DV.
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010.
%             http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf

%   MATPOWER
%   Copyright (c) 2008-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define
nl = length(mu);

diagmu = sparse(1:nl, 1:nl, mu, nl, nl);

[F11, F12, F21, F22] = d2F_dV2(V, conj(F) .* mu);
H11 = 2 * real( F11 + dF_dV1.' * diagmu * conj(dF_dV1) );
H21 = 2 * real( F21 + dF_dV2.' * diagmu * conj(dF_dV1) );
H12 = 2 * real( F12 + dF_dV1.' * diagmu * conj(dF_dV2) );
H22 = 2 * real( F22 + dF_dV2.' * diagmu * conj(dF_dV2) );
