function [Haa, Hav, Hva, Hvv] = d2Ibr_dV2(Ybr, V, lam)
%D2IBR_DV2   Computes 2nd derivatives of complex branch current w.r.t. voltage.
%   [HAA, HAV, HVA, HVV] = D2IBR_DV2(CBR, YBR, V, LAM) returns 4 matrices
%   containing the partial derivatives w.r.t. voltage angle and magnitude
%   of the product of a vector LAM with the 1st partial derivatives of the
%   complex branch currents. Takes sparse branch admittance matrix YBR,
%   voltage vector V and nl x 1 vector of multipliers LAM. Output matrices
%   are sparse.
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       Ybr = Yf;
%       [Haa, Hav, Hva, Hvv] = d2Ibr_dV2(Ybr, V, lam);
%
%   Here the output matrices correspond to:
%       Haa = (d/dVa (dIbr_dVa.')) * lam
%       Hav = (d/dVm (dIbr_dVa.')) * lam
%       Hva = (d/dVa (dIbr_dVm.')) * lam
%       Hvv = (d/dVm (dIbr_dVm.')) * lam
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
nb = length(V);

diaginvVm = sparse(1:nb, 1:nb, ones(nb, 1)./abs(V), nb, nb);

Haa = sparse(1:nb, 1:nb, -(Ybr.' * lam) .* V, nb, nb);
Hva = -1j * Haa * diaginvVm;
Hav = Hva;
Hvv = sparse(nb, nb);
