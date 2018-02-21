function [dSf_dVi, dSf_dVr, dSt_dVi, dSt_dVr, Sf, St] = dSbr_dV_C(branch, Yf, Yt, V)
%dSbr_dV_C   Computes partial derivatives of power flows w.r.t. voltage.
%   [DSF_DVI, DSF_DVR, DST_DVI, DST_DVR, SF, ST] = dSbr_dV_C(BRANCH, YF, YT, V)
%   returns four matrices containing partial derivatives of the complex
%   branch power flows at "from" and "to" ends of each branch w.r.t real
%   an imaginary part of complex voltage respectively (for all buses). If YF is a
%   sparse matrix, the partial derivative matrices will be as well. Optionally
%   returns vectors containing the power flows themselves. The following
%   explains the expressions used to form the matrices:
%
%   If = Yf * V;
%   Sf = diag(Vf) * conj(If) = diag(conj(If)) * Vf
%
%   Partials of V, Vf & If w.r.t. imaginary part of complex voltage
%       dV/dVi  = j * diag(ones(n,1))
%       dVf/dVi = j * Cf
%       dIf/dVi = Yf * j 
%   where Cf is connection matrix for line & from buses
%
%   Partials of V, Vf & If w.r.t. real part of complex voltage
%       dV/dVr  = diag(ones(n,1))
%       dVf/dVr = Cf
%       dIf/dVr = Yf 
%
%   Partials of Sf w.r.t. imaginary part of complex voltage
%       dSf/dVi = j * (conj(diag(If)) * Cf - diag(Vf) * conj(Yf))
%
%   Partials of Sf w.r.t. real part of complex voltage
%       dSf/dVr = conj(diag(If)) * Cf + diag(Vf) * conj(Yf)
%
%   Derivations for "to" bus are similar.
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dSf_dVi, dSf_dVr, dSt_dVi, dSt_dVr, Sf, St] = ...
%           dSbr_dV_C(branch, Yf, Yt, V);
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:

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

Cf = sparse(1:nl, f, ones(nl, 1), nl, nb);      %% connection matrix for line & from buses
Ct = sparse(1:nl, t, ones(nl, 1), nl, nb);      %% connection matrix for line & to buses    
if issparse(Yf)             %% sparse version (if Yf is sparse)
    diagVf      = sparse(1:nl, 1:nl, V(f), nl, nl);
    diagVt      = sparse(1:nl, 1:nl, V(t), nl, nl);
    diagIf      = sparse(1:nl, 1:nl, If, nl, nl);
    diagIt      = sparse(1:nl, 1:nl, It, nl, nl);
else                        %% dense version
    diagVf      = diag(V(f));
    diagVt      = diag(V(t));
    diagIf      = diag(If);    
    diagIt      = diag(It);
end
    dSf_dVi = 1j * (conj(diagIf)*Cf - diagVf * conj(Yf));   %% dSf_dVi
    dSf_dVr = conj(diagIf)*Cf + diagVf * conj(Yf) ;         %% dSf_dVr

    dSt_dVi = 1j * (conj(diagIt)*Ct - diagVt * conj(Yt));   %% dSt_dVi
    dSt_dVr = conj(diagIt)*Ct + diagVt * conj(Yt) ;         %% dSt_dVr
   
if nargout > 4
    Sf = V(f) .* conj(If);
    St = V(t) .* conj(It);
end
end
