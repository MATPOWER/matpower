function [G_SV, G_VS] = d2Imis_dVdSg_C(V, lam, Cg)
%d2Imis_dVdSg_C   Computes 2nd derivatives of current balance w.r.t. voltage and complex power generation.
%   [G_SV, G_VS] = d2Imis_dVdSg_C(V, LAM, CG) returns 2 matrices
%   containing the partial derivatives w.r.t. real and imaginary part of compelx volate, and real and reactive 
%   power generations of the product of a vector LAM with the 1st partial derivatives of the
%   complex current balance. Output matrices are sparse.
%
%   Example:
%       [G_SV, G_VS] = d2Imis_dVdSg_C(V, lam, Cg);
%
%   Here the output matrices correspond to:
%   G_SV = [ G_Pg_Vi G_Pg_Vr; 
%            G_Qg_Vi G_Qg_Vr ]
%   G_VS = G_SV.';
%   where
%       G_Pg_Vi = (d/dVi (dImis_dPg.')) * lam
%       G_Pg_Vr = (d/dVr (dImis_dPg.')) * lam
%       G_Qg_Vi = (d/dVi (dImis_dQg.')) * lam
%       G_Qg_Vr = (d/dVr (dImis_dQg.')) * lam
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER, see:

nb = length(V);
LamV2 = sparse(1:nb,1:nb, lam./(V.^2), nb, nb);
CgLamV2 = Cg'*conj(LamV2);
    
G_Pg_Vi = -1j*CgLamV2;
G_Pg_Vr = CgLamV2;
G_Qg_Vi = -CgLamV2;
G_Qg_Vr = -1j*CgLamV2;

G_SV = [ G_Pg_Vi G_Pg_Vr; 
         G_Qg_Vi G_Qg_Vr ];    
       
G_VS = G_SV.'; 
end