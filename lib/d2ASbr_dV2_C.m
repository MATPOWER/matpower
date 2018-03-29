function [Hrr, Hri, Hir, Hii] = ...
    d2ASbr_dV2_C(dSbr_dVr, dSbr_dVi, Sbr, Cbr, Ybr, V, mu)
%d2ASbr_dV2_C   Computes 2nd derivatives of |complex power flow|^2 w.r.t. V.
%   [HRR, HRI, HIR, HII] = d2ASbr_dV2_C(DSBR_DVR, DSBR_DVI, SBR, CBR, YBR, V, MU)
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
%       [dSf_dVr, dSf_dVi, dSt_dVr, dSt_dVi, Sf, St] = ...
%               dSbr_dV_C(branch, Yf, Yt, V);
%       Cbr = Cf;
%       Ybr = Yf;
%       dSbr_dVr = dSf_dVr;
%       dSbr_dVi = dSf_dVi;
%       Sbr = Sf;
%       [Hrr, Hri, Hir, Hii] = ...
%             d2ASbr_dV2_C(dSbr_dVr, dSbr_dVi, Sbr, Cbr, Ybr, V, mu);
%
%   Here the output matrices correspond to:
%     Hrr = (d/dVr (dASbr_dVr.')) * mu
%     Hri = (d/dVi (dASbr_dVr.')) * mu
%     Hir = (d/dVr (dASbr_dVi.')) * mu
%     Hii = (d/dVi (dASbr_dVi.')) * mu
%
%   See also DSBR_DV_C.
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010.
%             http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf

%   MATPOWER
%   Copyright (c) 2008-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define
nl = length(mu);

diagmu = sparse(1:nl, 1:nl, mu, nl, nl);
diagSbr_conj = sparse(1:nl, 1:nl, conj(Sbr), nl, nl);

[Srr, Sri, Sir, Sii] = d2Sbr_dV2_C(Cbr, Ybr, V, diagSbr_conj * mu);
Hrr = 2 * real( Srr + dSbr_dVr.' * diagmu * conj(dSbr_dVr) );
Hir = 2 * real( Sir + dSbr_dVi.' * diagmu * conj(dSbr_dVr) );
Hri = 2 * real( Sri + dSbr_dVr.' * diagmu * conj(dSbr_dVi) );
Hii = 2 * real( Sii + dSbr_dVi.' * diagmu * conj(dSbr_dVi) );
