function [H11, H12, H21, H22] = ...
    d2ASbr_dV2(dSbr_dV1, dSbr_dV2, Sbr, Cbr, Ybr, V, mu, vcart)
%D2ASBR_DV2   Computes 2nd derivatives of |power flow|^2 w.r.t. V.
%
%   -----  DEPRECATED - Please use D2ABR_DV2 instead    -----
%   -----  See wrapper code in D2ASBR_DV2 for example.  -----
%%
%   The derivatives can be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 8th argument.
%
%   [HAA, HAV, HVA, HVV] = D2ASBR_DV2(DSBR_DV1, DSBR_DV2, SBR, CBR, YBR, V, MU)
%   [HAA, HAV, HVA, HVV] = D2ASBR_DV2(DSBR_DV1, DSBR_DV2, SBR, CBR, YBR, V, MU, 0)
%
%   Returns 4 matrices containing the partial derivatives w.r.t. voltage
%   angle and magnitude of the product of a vector MU with the 1st partial
%   derivatives of the square of the magnitude of branch power flows.
%
%   [HRR, HRI, HIR, HII] = D2ASBR_DV2(DSBR_DV1, DSBR_DV2, SBR, CBR, YBR, V, MU, 1)
%
%   Returns 4 matrices containing the partial derivatives w.r.t. real and
%   imaginary part of complex voltage of the product of a vector MU with the
%   1st partial derivatives of the square of the magnitude of branch power
%   flows.
%
%   Takes as inputs sparse first derivative matrices of complex flow, complex
%   flow vector, sparse connection matrix CBR, sparse branch admittance matrix
%   YBR, voltage vector V and nl x 1 vector of multipliers MU. Output matrices
%   are sparse.
%
%   Example:
%       f = branch(:, F_BUS);
%       Cf =  sparse(1:nl, f, ones(nl, 1), nl, nb);
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = ...
%               dSbr_dV(branch, Yf, Yt, V, vcart);
%       Cbr = Cf;
%       Ybr = Yf;
%       dSbr_dV1 = dSf_dV1;
%       dSbr_dV2 = dSf_dV2;
%       Sbr = Sf;
%       [H11, H12, H21, H22] = ...
%             d2ASbr_dV2(dSbr_dV1, dSbr_dV2, Sbr, Cbr, Ybr, V, mu, vcart);
%
%   Here the output matrices correspond to:
%     H11 = d/dV1 (dASbr_dV1.' * mu)
%     H12 = d/dV2 (dASbr_dV1.' * mu)
%     H21 = d/dV1 (dASbr_dV2.' * mu)
%     H22 = d/dV2 (dASbr_dV2.' * mu)
%
%   See also DSBR_DV, DABR_DV.
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
if nargin < 8
    vcart = 0;      %% default to polar coordinates
end

d2F_dV2 = @(V, mu)d2Sbr_dV2(Cbr, Ybr, V, mu, vcart);
[H11, H12, H21, H22] = d2Abr_dV2(d2F_dV2, dSbr_dV1, dSbr_dV2, Sbr, V, mu);
