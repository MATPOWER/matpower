function [H11, H12, H21, H22] = ...
    d2AIbr_dV2(dIbr_dV1, dIbr_dV2, Ibr, Ybr, V, mu, vcart)
%D2AIBR_DV2   Computes 2nd derivatives of |complex current|^2 w.r.t. V.
%
%   -----  DEPRECATED - Please use D2ABR_DV2 instead    -----
%   -----  See wrapper code in D2AIBR_DV2 for example.  -----
%
%   The derivatives can be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 8th argument.
%
%   [HAA, HAV, HVA, HVV] = D2AIBR_DV2(DIBR_DVA, DIBR_DVM, IBR, YBR, V, MU)
%   [HAA, HAV, HVA, HVV] = D2AIBR_DV2(DIBR_DVA, DIBR_DVM, IBR, YBR, V, MU, 0)
%
%   Returns 4 matrices containing the partial derivatives w.r.t. voltage
%   angle and magnitude of the product of a vector MU with the 1st partial
%   derivatives of the square of the magnitude the branch currents.
%
%   [HRR, HRI, HIR, HII] = D2AIBR_DV2(DIBR_DVA, DIBR_DVM, IBR, YBR, V, MU, 1)
%
%   Returns 4 matrices containing the partial derivatives w.r.t. real and
%   imaginary part of complex voltage of the product of a vector MU with the
%   1st partial derivatives of the square of the magnitude the branch currents.
%
%   Takes as inputs sparse first derivative matrices of complex current,
%   complex current vector, sparse branch admittance matrix YBR, voltage
%   vector V and nl x 1 vector of multipliers MU. Output matrices are sparse.
%
%   Example:
%       f = branch(:, F_BUS);
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dIf_dV1, dIf_dV2, dIt_dV1, dIt_dV2, If, It] = ...
%               dIbr_dV(branch, Yf, Yt, V, vcart);
%       Ybr = Yf;
%       dIbr_dV1 = dIf_dV1;
%       dIbr_dV2 = dIf_dV2;
%       Ibr = If;
%       [H11, H12, H21, H22] = ...
%             d2AIbr_dV2(dIbr_dV1, dIbr_dV2, Ibr, Ybr, V, mu, vcart);
%
%   Here the output matrices correspond to:
%     H11 = d/dV1 (dAIbr_dV1.' * mu)
%     H12 = d/dV2 (dAIbr_dV1.' * mu)
%     H21 = d/dV1 (dAIbr_dV2.' * mu)
%     H22 = d/dV2 (dAIbr_dV2.' * mu)
%
%   See also DIBR_DV, DABR_DV.
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

%% default input args
if nargin < 7
    vcart = 0;      %% default to polar coordinates
end

d2F_dV2 = @(V, mu)d2Ibr_dV2(Ybr, V, mu, vcart);
[H11, H12, H21, H22] = d2Abr_dV2(d2F_dV2, dIbr_dV1, dIbr_dV2, Ibr, V, mu);
