function [Hii, Hir, Hri, Hrr] = ...
    d2ASbr_dV2_C(dSbr_dVi, dSbr_dVr, Sbr, Cbr, Ybr, V, mu)
%d2ASbr_dV2_C   Computes 2nd derivatives of |complex power flow|^2 w.r.t. V.
%   [HII, HIR, HRI, HRR] = d2ASbr_dV2_C(DSBR_DVI, DSBR_DVR, SBR, CBR, YBR, V, MU)
%   returns 4 matrices containing the partial derivatives w.r.t. real and
%   imaginary part of complex voltage of the product of a vector MU with the 1st partial
%   derivatives of the square of the magnitude of branch complex power flows.
%   Takes sparse first derivative matrices of complex flow, complex flow
%   vector, sparse connection matrix CBR, sparse branch admittance matrix YBR,
%   voltage vector V and nl x 1 vector of multipliers MU. Output matrices
%   are sparse.
%
%   Example:
%       f = branch(:, F_BUS);
%       Cf =  sparse(1:nl, f, ones(nl, 1), nl, nb);
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dSf_dVi, dSf_dVr, dSt_dVi, dSt_dVr, Sf, St] = ...
%               dSbr_dV_C(branch, Yf, Yt, V);
%       Cbr = Cf;
%       Ybr = Yf;
%       dSbr_dVi = dSf_dVi;
%       dSbr_dVr = dSf_dVr;
%       Sbr = Sf;
%       [Hii, Hir, Hri, Hrr] = ...
%             d2ASbr_dV2_C(dSbr_dVi, dSbr_dVr, Sbr, Cbr, Ybr, V, mu);
%
%   Here the output matrices correspond to:
%     Hii = (d/dVi (dASbr_dVi.')) * mu
%     Hir = (d/dVr (dASbr_dVi.')) * mu
%     Hri = (d/dVi (dASbr_dVr.')) * mu
%     Hrr = (d/dVr (dASbr_dVr.')) * mu
%
%   See also DSBR_DV_C.
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:

%% define
nl = length(mu);

diagmu = sparse(1:nl, 1:nl, mu, nl, nl);
diagSbr_conj = sparse(1:nl, 1:nl, conj(Sbr), nl, nl);

[Sii, Sir, Sri, Srr] = d2Sbr_dV2_C(Cbr, Ybr, V, diagSbr_conj * mu);
Hii = 2 * real( Sii + dSbr_dVi.' * diagmu * conj(dSbr_dVi) );
Hri = 2 * real( Sri + dSbr_dVr.' * diagmu * conj(dSbr_dVi) );
Hir = 2 * real( Sir + dSbr_dVi.' * diagmu * conj(dSbr_dVr) );
Hrr = 2 * real( Srr + dSbr_dVr.' * diagmu * conj(dSbr_dVr) );
