function [Hii, Hir, Hri, Hrr] = d2Sbr_dV2_C(Cbr, Ybr, V, mu)
%d2Sbr_dV2_C   Computes 2nd derivatives of complex power flow w.r.t. voltage.
%   [HII, HIR, HRI, HRR] = d2Sbr_dV2_C(CBR, YBR, V, MU) returns 4 matrices
%   containing the partial derivatives w.r.t. real and imaginary part of
%   complex voltage of the product of a vector MU with the 1st partial derivatives of the
%   complex branch power flows. Takes sparse connection matrix CBR, sparse
%   branch admittance matrix YBR, voltage vector V and nl x 1 vector of
%   multipliers MU. Output matrices are sparse.
%
%   Example:
%       f = branch(:, F_BUS);
%       Cf =  sparse(1:nl, f, ones(nl, 1), nl, nb);
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       Cbr = Cf;
%       Ybr = Yf;
%       [Hii, Hir, Hri, Hrr] = d2Sbr_dV2_C(Cbr, Ybr, V, mu);
%
%   Here the output matrices correspond to:
%       Hii = (d/dVi (dSbr_dVi.')) * mu
%       Hir = (d/dVr (dSbr_dVi.')) * mu
%       Hri = (d/dVi (dSbr_dVr.')) * mu
%       Hrr = (d/dVr (dSbr_dVr.')) * mu
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:

%% define
nl = length(mu);
diagSMu = sparse(1:nl, 1:nl, mu, nl, nl);

CmuY    = Cbr'*diagSMu*conj(Ybr);
YmuC    = Ybr'*diagSMu*Cbr;

Hii = CmuY + YmuC;
Hir = 1j*(CmuY - YmuC);
Hrr = Hii;
Hri = -Hir;    
end
