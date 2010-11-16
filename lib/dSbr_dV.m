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
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2010 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

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
