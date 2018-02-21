function [dAf_dVi, dAf_dVr, dAt_dVi, dAt_dVr] = ...
                        dAbr_dV_C(dSf_dVi, dSf_dVr, dSt_dVi, dSt_dVr, Sf, St)
%dAbr_dV_C   Partial derivatives of squared flow magnitudes w.r.t voltage.
%   [DAF_DVI, DAF_DVR, DAT_DVI, DAT_DVR] = ...
%               dAbr_dV_C(DFF_DVI, DFF_DVR, DFT_DVI, DFT_DVR, FF, FT)
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
%   Partial w.r.t Vm & Va
%       dAf/dVr = dAf/dPf * dPf/dVr + dAf/dQf * dQf/dVr
%       dAf/dVi = dAf/dPf * dPf/dVi + dAf/dQf * dQf/dVi
%
%   Derivations for "to" bus are similar.
%
%   Examples:
%       %% squared apparent power flow
%       [dFf_dVi, dFf_dVr, dFt_dVi, dFt_dVr, Ff, Ft] = ...
%               dSbr_dV_C(branch(il,:), Yf, Yt, V);
%       [dAf_dVi, dAf_dVr, dAt_dVi, dAt_dVr] = ...
%               dAbr_dV_C(dFf_dVi, dFf_dVr, dFt_dVi, dFt_dVr, Ff, Ft);
%
%   See also DIBR_DV_C, DSBR_DV_C.
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:
                  
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
end
