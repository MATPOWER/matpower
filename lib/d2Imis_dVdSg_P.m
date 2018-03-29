function [G_SV, G_VS] = d2Imis_dVdSg_P(V, lam, Cg)
%d2Imis_dVdSg_P   Computes 2nd derivatives of current balance w.r.t. voltage and complex power generation.
%   [G_SV, G_VS] = d2Imis_dVdSg_P(V, LAM, CG) returns 2 matrices
%   containing the partial derivatives w.r.t. voltage magnitude, voltage angle, and real and reactive
%   power generations of the product of a vector LAM with the 1st partial derivatives of the
%   complex current balance. Output matrices are sparse.
%
%   Example:
%       [G_SV, G_VS] = d2Imis_dVdSg_P(V, lam, Cg);
%
%   Here the output matrices correspond to:
%   G_SV = [ G_Pg_Va G_Pg_Vm;
%            G_Qg_Va G_Qg_Vm ]
%   G_VS = G_SV.';
%   where
%       G_Pg_Va = (d/dVa (dImis_dPg.')) * lam
%       G_Pg_Vm = (d/dVm (dImis_dPg.')) * lam
%       G_Qg_Va = (d/dVa (dImis_dQg.')) * lam
%       G_Qg_Vm = (d/dVm (dImis_dQg.')) * lam
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER, see:

nb = length(V);
diagVm1     = sparse(1:nb, 1:nb, 1./abs(V), nb, nb);
LamConjV1   = sparse(1:nb, 1:nb, lam./conj(V), nb, nb);
CgLamConjV1 = Cg'*LamConjV1;

G_Pg_Va = -1j*CgLamConjV1;
G_Pg_Vm = CgLamConjV1*diagVm1;
G_Qg_Va = -CgLamConjV1;
G_Qg_Vm = -1j*CgLamConjV1*diagVm1;

G_SV = [ G_Pg_Va G_Pg_Vm;
         G_Qg_Va G_Qg_Vm ];

G_VS = G_SV.';
