function [dP_dV, dP_dlam] = cpf_p_jac(parameterization, z, V, lam, Vprv, lamprv, pv, pq)
%CPF_P_JAC Computes partial derivatives of CPF parameterization function.
%   [DP_DV, DP_DLAM ] = CPF_P_JAC(PARAMETERIZATION, Z, V, LAM, ...
%                                                   VPRV, LAMPRV, PV, PQ)
%
%   Computes the partial derivatives of the continuation power flow
%   parameterization function w.r.t. bus voltages and the continuation
%   parameter lambda.
%
%   Inputs:
%       PARAMETERIZATION : Value of cpf.parameterization option.
%       Z : normalized tangent prediction vector from previous step
%       V : complex bus voltage vector at current solution
%       LAM : scalar lambda value at current solution
%       VPRV : complex bus voltage vector at previous solution
%       LAMPRV : scalar lambda value at previous solution
%       PV : vector of indices of PV buses
%       PQ : vector of indices of PQ buses
%
%   Outputs:
%       DP_DV : partial of parameterization function w.r.t. voltages
%       DP_DLAM : partial of parameterization function w.r.t. lambda
%
%   See also CPF_PREDICTOR, CPF_CORRECTOR.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Shrirang Abhyankar, Argonne National Laboratory
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if parameterization == 1        %% natural
    npv = length(pv);
    npq = length(pq);
    dP_dV = zeros(1, npv+2*npq);
    if lam >= lamprv
        dP_dlam = 1.0;
    else
        dP_dlam = -1.0;
    end
elseif parameterization == 2    %% arc length
    Va = angle(V);
    Vm = abs(V);
    Vaprv = angle(Vprv);
    Vmprv = abs(Vprv);
    dP_dV = 2*([Va([pv; pq]); Vm(pq)] - [Vaprv([pv; pq]); Vmprv(pq)])';
    if lam == lamprv    %% first step
        dP_dlam = 1.0;  %% avoid singular Jacobian that would result
                        %% from [dP_dV, dP_dlam] = 0
    else
        dP_dlam = 2*(lam-lamprv);
    end
elseif parameterization == 3    %% pseudo arc length
    nb = length(V);
    dP_dV = z([pv; pq; nb+pq])';
    dP_dlam = z(2*nb+1);
end 
