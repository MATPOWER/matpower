function P = cpf_p(parameterization, step, z, V, lam, Vprv, lamprv, pv, pq)
%CPF_P Computes the value of the CPF parameterization function.
%   P = CPF_P(PARAMETERIZATION, STEP, Z, V, LAM, VPRV, LAMPRV, PV, PQ)
%
%   Computes the value of the parameterization function at the current
%   solution point.
%
%   Inputs:
%       PARAMETERIZATION : Value of cpf.parameterization option
%       STEP : continuation step size
%       Z : normalized tangent prediction vector from previous step
%       V : complex bus voltage vector at current solution
%       LAM : scalar lambda value at current solution
%       VPRV : complex bus voltage vector at previous solution
%       LAMPRV : scalar lambda value at previous solution
%       PV : vector of indices of PV buses
%       PQ : vector of indices of PQ buses
%
%   Outputs:
%       P : value of the parameterization function at the current point
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

%% evaluate P(x0, lambda0)
if parameterization == 1        %% natural
    if lam >= lamprv
        P = lam - lamprv - step;
    else
        P = lamprv - lam - step;
    end
elseif parameterization == 2    %% arc length
    Va = angle(V);
    Vm = abs(V);
    Vaprv = angle(Vprv);
    Vmprv = abs(Vprv);
    P = sum(([Va([pv; pq]); Vm(pq); lam] - [Vaprv([pv; pq]); Vmprv(pq); lamprv]).^2) - step^2;
elseif parameterization == 3    %% pseudo arc length
    nb = length(V);
    Va = angle(V);
    Vm = abs(V);
    Vaprv = angle(Vprv);
    Vmprv = abs(Vprv);
    P = z([pv; pq; nb+pq; 2*nb+1])' * ...
        ( [Va([pv; pq]); Vm(pq); lam] - [Vaprv([pv; pq]); Vmprv(pq); lamprv] )...
        - step;
end 
