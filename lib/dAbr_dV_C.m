function [dAf_dVr, dAf_dVi, dAt_dVr, dAt_dVi] = ...
                        dAbr_dV_C(dSf_dVr, dSf_dVi, dSt_dVr, dSt_dVi, Sf, St)
%dAbr_dV_C   Partial derivatives of squared flow magnitudes w.r.t voltage.
%   [DAF_DVR, DAF_DVI, DAT_DVR, DAT_DVI] = ...
%               dAbr_dV_C(DFF_DVR, DFF_DVI, DFT_DVR, DFT_DVI, FF, FT)
%   returns four matrices containing partial derivatives of the square of
%   the branch flow magnitudes at "from" & "to" ends of each branch w.r.t
%   real and imaginary part of complex voltage respectively (for all buses), given
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
%   Partial w.r.t Vr & Vi
%       dAf/dVr = dAf/dPf * dPf/dVr + dAf/dQf * dQf/dVr
%       dAf/dVi = dAf/dPf * dPf/dVi + dAf/dQf * dQf/dVi
%
%   Derivations for "to" bus are similar.
%
%   Examples:
%       %% squared apparent power flow
%       [dFf_dVr, dFf_dVi, dFt_dVr, dFt_dVi, Ff, Ft] = ...
%               dSbr_dV_C(branch(il,:), Yf, Yt, V);
%       [dAf_dVr, dAf_dVi, dAt_dVr, dAt_dVi] = ...
%               dAbr_dV_C(dFf_dVr, dFf_dVi, dFt_dVr, dFt_dVi, Ff, Ft);
%
%   See also DIBR_DV_C, DSBR_DV_C.
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010.
%             http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% dimensions
nl = length(Sf);

%%----- partials w.r.t. real and reactive power flows -----
dAf_dPf = sparse(1:nl, 1:nl, 2 * real(Sf), nl, nl);
dAf_dQf = sparse(1:nl, 1:nl, 2 * imag(Sf), nl, nl);
dAt_dPt = sparse(1:nl, 1:nl, 2 * real(St), nl, nl);
dAt_dQt = sparse(1:nl, 1:nl, 2 * imag(St), nl, nl);

%% partials w.r.t. voltage magnitudes and angles
dAf_dVr = dAf_dPf * real(dSf_dVr) + dAf_dQf * imag(dSf_dVr);
dAf_dVi = dAf_dPf * real(dSf_dVi) + dAf_dQf * imag(dSf_dVi);
dAt_dVr = dAt_dPt * real(dSt_dVr) + dAt_dQt * imag(dSt_dVr);
dAt_dVi = dAt_dPt * real(dSt_dVi) + dAt_dQt * imag(dSt_dVi);
