function [dAf_dV1, dAf_dV2, dAt_dV1, dAt_dV2] = ...
                        dAbr_dV(dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2, Ff, Ft)
%DABR_DV   Partial derivatives of squared flow magnitudes w.r.t voltage.
%   [DAF_DV1, DAF_DV2, DAT_DV1, DAT_DV2] = ...
%               DABR_DV(DFF_DV1, DFF_DV2, DFT_DV1, DFT_DV2, FF, FT)
%   returns four matrices containing partial derivatives of the square of
%   the branch flow magnitudes at "from" & "to" ends of each branch w.r.t
%   voltage components (either angle and magnitude, respectively, if polar,
%   or real and imaginary, respectively, if cartesian) for all buses, given
%   the flows and flow sensitivities. Flows could be complex current or
%   complex or real power. Notation below is based on complex power. The
%   following explains the expressions used to form the matrices:
%
%   Let Af refer to the square of the apparent power at the "from" end of
%   each branch,
%
%       Af = abs(Sf).^2
%          = Sf .* conj(Sf)
%          = Pf.^2 + Qf.^2
%
%   then ...
%
%   Partial w.r.t real power,
%       dAf/dPf = 2 * diag(Pf)
%
%   Partial w.r.t reactive power,
%       dAf/dQf = 2 * diag(Qf)
%
%   Partial w.r.t V1 & V2 (e.g. Va and Vm, or Vr and Vi)
%       dAf/dV1 = dAf/dPf * dPf/dV1 + dAf/dQf * dQf/dV1
%       dAf/dV2 = dAf/dPf * dPf/dV2 + dAf/dQf * dQf/dV2
%
%   Derivations for "to" bus are similar.
%
%   Examples:
%       %% squared current magnitude
%       [dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2, Ff, Ft] = ...
%               dIbr_dV(branch(il,:), Yf, Yt, V);
%       [dAf_dV1, dAf_dV2, dAt_dV1, dAt_dV2] = ...
%               dAbr_dV(dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2, Ff, Ft);
%
%       %% squared apparent power flow
%       [dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2, Ff, Ft] = ...
%               dSbr_dV(branch(il,:), Yf, Yt, V);
%       [dAf_dV1, dAf_dV2, dAt_dV1, dAt_dV2] = ...
%               dAbr_dV(dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2, Ff, Ft);
%
%       %% squared real power flow
%       [dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2, Ff, Ft] = ...
%               dSbr_dV(branch(il,:), Yf, Yt, V);
%       dFf_dV1 = real(dFf_dV1);
%       dFf_dV2 = real(dFf_dV2);
%       dFt_dV1 = real(dFt_dV1);
%       dFt_dV2 = real(dFt_dV2);
%       [dAf_dV1, dAf_dV2, dAt_dV1, dAt_dV2] = ...
%               dAbr_dV(dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2, Ff, Ft);
%
%   See also DIBR_DV, DSBR_DV.
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010.
%             http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf

%   MATPOWER
%   Copyright (c) 1996-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% dimensions
nl = length(Ff);

%%----- partials w.r.t. real and imaginary flows -----
dAf_dFfr = sparse(1:nl, 1:nl, 2 * real(Ff), nl, nl);
dAf_dFfi = sparse(1:nl, 1:nl, 2 * imag(Ff), nl, nl);
dAt_dFtr = sparse(1:nl, 1:nl, 2 * real(Ft), nl, nl);
dAt_dFti = sparse(1:nl, 1:nl, 2 * imag(Ft), nl, nl);

%% partials w.r.t. voltage components (angle, magnitude or real, imaginary)
dAf_dV1 = dAf_dFfr * real(dFf_dV1) + dAf_dFfi * imag(dFf_dV1);
dAf_dV2 = dAf_dFfr * real(dFf_dV2) + dAf_dFfi * imag(dFf_dV2);
dAt_dV1 = dAt_dFtr * real(dFt_dV1) + dAt_dFti * imag(dFt_dV1);
dAt_dV2 = dAt_dFtr * real(dFt_dV2) + dAt_dFti * imag(dFt_dV2);
