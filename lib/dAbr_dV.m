function [dAf_dVa, dAf_dVm, dAt_dVa, dAt_dVm] = ...
                        dAbr_dV(dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St)
%DABR_DV   Partial derivatives of squared flow magnitudes w.r.t voltage.
%   [DAF_DVA, DAF_DVM, DAT_DVA, DAT_DVM] = ...
%               DABR_DV(DFF_DVA, DFF_DVM, DFT_DVA, DFT_DVM, FF, FT)
%   returns four matrices containing partial derivatives of the square of
%   the branch flow magnitudes at "from" & "to" ends of each branch w.r.t
%   voltage magnitude and voltage angle respectively (for all buses), given
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
%   Partial w.r.t Vm & Va
%       dAf/dVm = dAf/dPf * dPf/dVm + dAf/dQf * dQf/dVm
%       dAf/dVa = dAf/dPf * dPf/dVa + dAf/dQf * dQf/dVa
%
%   Derivations for "to" bus are similar.
%
%   Examples:
%       %% squared current magnitude
%       [dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft] = ...
%               dIbr_dV(branch(il,:), Yf, Yt, V);
%       [dAf_dVa, dAf_dVm, dAt_dVa, dAt_dVm] = ...
%               dAbr_dV(dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft);
%
%       %% squared apparent power flow
%       [dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft] = ...
%               dSbr_dV(branch(il,:), Yf, Yt, V);
%       [dAf_dVa, dAf_dVm, dAt_dVa, dAt_dVm] = ...
%               dAbr_dV(dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft);
%
%       %% squared real power flow
%       [dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft] = ...
%               dSbr_dV(branch(il,:), Yf, Yt, V);
%       dFf_dVa = real(dFf_dVa);
%       dFf_dVm = real(dFf_dVm);
%       dFt_dVa = real(dFt_dVa);
%       dFt_dVm = real(dFt_dVm);
%       [dAf_dVa, dAf_dVm, dAt_dVa, dAt_dVm] = ...
%               dAbr_dV(dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft);
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
dAf_dVm = dAf_dPf * real(dSf_dVm) + dAf_dQf * imag(dSf_dVm);
dAf_dVa = dAf_dPf * real(dSf_dVa) + dAf_dQf * imag(dSf_dVa);
dAt_dVm = dAt_dPt * real(dSt_dVm) + dAt_dQt * imag(dSt_dVm);
dAt_dVa = dAt_dPt * real(dSt_dVa) + dAt_dQt * imag(dSt_dVa);
