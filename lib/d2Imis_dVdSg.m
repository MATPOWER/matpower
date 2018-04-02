function Gsv = d2Imis_dVdSg(Cg, V, lam, vcart)
%D2IMIS_DVDSG   Computes 2nd derivatives of current balance w.r.t. V and Sg.
%
%   The derivatives can be take with respect to polar or cartesian coordinates
%   of voltage, depending on the 4th argument.
%
%   GSV = D2IMIS_DVDSG(CG, V, LAM)
%   GSV = D2IMIS_DVDSG(CG, V, LAM, 0)
%
%   Returns a matrix containing the partial derivatives w.r.t. voltage angle
%   and magnitude of the product of a vector LAM with the 1st partial
%   derivatives of the real and reactive power generation.
%
%   GSV = D2IMIS_DVDSG(CG, V, LAM, 1)
%
%   Returns a matrix containing the partial derivatives w.r.t. real and
%   imaginary parts of voltage of the product of a vector LAM with the 1st
%   partial derivatives of the real and reactive power generation.
%
%   Takes the generator connection matrix, complex voltage vector V and
%   nb x 1 vector of multipliers LAM. Output matrices are sparse.
%
%   Examples:
%       Cg = sparse(gen(:, GEN_BUS), 1:ng, -, nb, ng);
%       Gsv = d2Imis_dVdSg(Cg, V, lam);
%
%       Here the output matrix corresponds to:
%           Gsv = [ Gpa Gpv;
%                   Gqa Gqv ];
%           Gpa = d/dVa (dImis_dPg.' * lam)
%           Gpv = d/dVm (dImis_dPg.' * lam)
%           Gqa = d/dVa (dImis_dQg.' * lam)
%           Gqv = d/dVm (dImis_dQg.' * lam)
%
%       [Grr, Gri, Gir, Gii] = d2Imis_dVdSg(Cg, V, lam, 1);
%
%       Here the output matrices correspond to:
%           Gsv = [ Gpr Gpi;
%                   Gqr Gqi ];
%           Gpr = d/dVr (dImis_dPg.' * lam)
%           Gpi = d/dVi (dImis_dPg.' * lam)
%           Gqr = d/dVr (dImis_dQg.' * lam)
%           Gqi = d/dVi (dImis_dQg.' * lam)
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010.
%             http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf
%   [TN3]  B. Sereeter and R. D. Zimmerman, "Addendum to AC Power Flows and
%          their Derivatives using Complex Matrix Notation: Nodal Current
%          Balance," MATPOWER Technical Note 3, April 2018.
%             http://www.pserc. cornell.edu/matpower/
%                                           TN3-More-OPF-Derivatives.pdf
%   [TN4]  B. Sereeter and R. D. Zimmerman, "AC Power Flows and their
%          Derivatives using Complex Matrix Notation and Cartesian
%          Coordinate Voltages," MATPOWER Technical Note 4, April 2018.
%             http://www.pserc.cornell.edu/matpower/
%                                           TN4-OPF-Derivatives-Cartesian.pdf

%   MATPOWER
%   Copyright (c) 2018, Power Systems Engineering Research Center (PSERC)
%   by Baljinnyam Sereeter, Delft University of Technology
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default input args
if nargin < 4
    vcart = 0;      %% default to polar coordinates
end

nb = length(V);

if vcart
    D = Cg' * sparse(1:nb,1:nb, lam./conj(V.^2), nb, nb);

    G_Pg_Vr = D;
    G_Pg_Vi = -1j * D;
    G_Qg_Vr = G_Pg_Vi;
    G_Qg_Vi = -D;

    Gsv = [ G_Pg_Vr G_Pg_Vi;
            G_Qg_Vr G_Qg_Vi ];
else
    D = sparse(1:nb, 1:nb, 1./abs(V), nb, nb);
    E = sparse(1:nb, 1:nb, lam./conj(V), nb, nb);
    K = Cg' * E;
    L = K * D;

    G_Pg_Va = -1j * K;
    G_Pg_Vm = L;
    G_Qg_Va = -K;
    G_Qg_Vm = -1j * L;

    Gsv = [ G_Pg_Va G_Pg_Vm;
            G_Qg_Va G_Qg_Vm ];
end
