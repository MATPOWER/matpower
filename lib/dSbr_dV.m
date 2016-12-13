function [dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = dSbr_dV(branch, Yf, Yt, V)
%DSBR_DV   Computes partial derivatives of power flows w.r.t. voltage.
%   [DSF_DVA, DSF_DVM, DST_DVA, DST_DVM, SF, ST] = DSBR_DV(BRANCH, YF, YT, V)
%   returns four matrices containing partial derivatives of the complex
%   branch power flows at "from" and "to" ends of each branch w.r.t voltage
%   magnitude and voltage angle respectively (for all buses). If YF is a
%   sparse matrix, the partial derivative matrices will be as well. Optionally
%   returns vectors containing the power flows themselves. The following
%   explains the expressions used to form the matrices:
%
%   If = Yf * V;
%   Sf = diag(Vf) * conj(If) = diag(conj(If)) * Vf
%
%   Partials of V, Vf & If w.r.t. voltage angles
%       dV/dVa  = j * diag(V)
%       dVf/dVa = sparse(1:nl, f, j * V(f)) = j * sparse(1:nl, f, V(f))
%       dIf/dVa = Yf * dV/dVa = Yf * j * diag(V)
%
%   Partials of V, Vf & If w.r.t. voltage magnitudes
%       dV/dVm  = diag(V./abs(V))
%       dVf/dVm = sparse(1:nl, f, V(f)./abs(V(f))
%       dIf/dVm = Yf * dV/dVm = Yf * diag(V./abs(V))
%
%   Partials of Sf w.r.t. voltage angles
%       dSf/dVa = diag(Vf) * conj(dIf/dVa)
%                       + diag(conj(If)) * dVf/dVa
%               = diag(Vf) * conj(Yf * j * diag(V))
%                       + conj(diag(If)) * j * sparse(1:nl, f, V(f))
%               = -j * diag(Vf) * conj(Yf * diag(V))
%                       + j * conj(diag(If)) * sparse(1:nl, f, V(f))
%               = j * (conj(diag(If)) * sparse(1:nl, f, V(f))
%                       - diag(Vf) * conj(Yf * diag(V)))
%
%   Partials of Sf w.r.t. voltage magnitudes
%       dSf/dVm = diag(Vf) * conj(dIf/dVm)
%                       + diag(conj(If)) * dVf/dVm
%               = diag(Vf) * conj(Yf * diag(V./abs(V)))
%                       + conj(diag(If)) * sparse(1:nl, f, V(f)./abs(V(f)))
%
%   Derivations for "to" bus are similar.
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = ...
%           dSbr_dV(branch, Yf, Yt, V);
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

%% define named indices into bus, gen, branch matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% define
f = branch(:, F_BUS);       %% list of "from" buses
t = branch(:, T_BUS);       %% list of "to" buses
nl = length(f);
nb = length(V);

%% compute currents
If = Yf * V;
It = Yt * V;

Vnorm = V ./ abs(V);
if issparse(Yf)             %% sparse version (if Yf is sparse)
    diagVf      = sparse(1:nl, 1:nl, V(f), nl, nl);
    diagIf      = sparse(1:nl, 1:nl, If, nl, nl);
    diagVt      = sparse(1:nl, 1:nl, V(t), nl, nl);
    diagIt      = sparse(1:nl, 1:nl, It, nl, nl);
    diagV       = sparse(1:nb, 1:nb, V, nb, nb);
    diagVnorm   = sparse(1:nb, 1:nb, Vnorm, nb, nb);
    
    dSf_dVa = 1j * (conj(diagIf) * sparse(1:nl, f, V(f), nl, nb) - diagVf * conj(Yf * diagV));
    dSf_dVm = diagVf * conj(Yf * diagVnorm) + conj(diagIf) * sparse(1:nl, f, Vnorm(f), nl, nb);
    dSt_dVa = 1j * (conj(diagIt) * sparse(1:nl, t, V(t), nl, nb) - diagVt * conj(Yt * diagV));
    dSt_dVm = diagVt * conj(Yt * diagVnorm) + conj(diagIt) * sparse(1:nl, t, Vnorm(t), nl, nb);
else                        %% dense version
    diagVf      = diag(V(f));
    diagIf      = diag(If);
    diagVt      = diag(V(t));
    diagIt      = diag(It);
    diagV       = diag(V);
    diagVnorm   = diag(Vnorm);
    temp1       = zeros(nl, nb);    temp1(sub2ind([nl,nb], (1:nl)', f)) = V(f);
    temp2       = zeros(nl, nb);    temp2(sub2ind([nl,nb], (1:nl)', f)) = Vnorm(f);
    temp3       = zeros(nl, nb);    temp3(sub2ind([nl,nb], (1:nl)', t)) = V(t);
    temp4       = zeros(nl, nb);    temp4(sub2ind([nl,nb], (1:nl)', t)) = Vnorm(t);
    
    dSf_dVa = 1j * (conj(diagIf) * temp1 - diagVf * conj(Yf * diagV));
    dSf_dVm = diagVf * conj(Yf * diagVnorm) + conj(diagIf) * temp2;
    dSt_dVa = 1j * (conj(diagIt) * temp3 - diagVt * conj(Yt * diagV));
    dSt_dVm = diagVt * conj(Yt * diagVnorm) + conj(diagIt) * temp4;
end

if nargout > 4
    Sf = V(f) .* conj(If);
    St = V(t) .* conj(It);
end
