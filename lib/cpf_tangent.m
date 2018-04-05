function z = cpf_tangent(V, lam, Ybus, Sbusb, Sbust, pv, pq, ...
                            zprv, Vprv, lamprv, parameterization, direction)
%CPF_TANGENT  Computes normalized tangent predictor for continuation power flow
%   Z = CPF_TANGENT(V, LAM, YBUS, SBUSB, SBUST, PV, PQ, ...
%                                 ZPRV, VPRV, LAMPRV, PARAMETERIZATION, DIRECTION)
%
%   Computes a normalized tangent predictor for the continuation power flow.
%
%   Inputs:
%       V : complex bus voltage vector at current solution
%       LAM : scalar lambda value at current solution
%       YBUS : complex bus admittance matrix
%       SBUSB : handle of function returning nb x 1 vector of complex
%               base case injections in p.u. and derivatives w.r.t. |V|
%       SBUST : handle of function returning nb x 1 vector of complex
%               target case injections in p.u. and derivatives w.r.t. |V|
%       PV : vector of indices of PV buses
%       PQ : vector of indices of PQ buses
%       ZPRV : normalized tangent prediction vector from previous step
%       VPRV : complex bus voltage vector at previous solution
%       LAMPRV : scalar lambda value at previous solution
%       PARAMETERIZATION : value of cpf.parameterization option.
%       DIRECTION: continuation direction (+1 for postive lambda
%                  increase, -1 otherwise)
%
%   Outputs:
%       Z : the normalized tangent prediction vector

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Shrirang Abhyankar, Argonne National Laboratory
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% sizes
nb = length(V);
npv = length(pv);
npq = length(pq);

Vm = abs(V);

%% compute Jacobian for the power flow equations
[dSbus_dVa, dSbus_dVm] = dSbus_dV(Ybus, V);
[dummy, neg_dSdb_dVm] = Sbusb(Vm);
[dummy, neg_dSdt_dVm] = Sbust(Vm);
dSbus_dVm = dSbus_dVm - neg_dSdb_dVm - lam * (neg_dSdt_dVm - neg_dSdb_dVm);

j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
j12 = real(dSbus_dVm([pv; pq], pq));
j21 = imag(dSbus_dVa(pq, [pv; pq]));
j22 = imag(dSbus_dVm(pq, pq));

J = [   j11 j12;
        j21 j22;    ];

Sxf = Sbust(Vm) - Sbusb(Vm);    %% "transfer" at current voltage level
dF_dlam = -[real(Sxf([pv; pq])); imag(Sxf(pq))];
[dP_dV, dP_dlam] = cpf_p_jac(parameterization, zprv, V, lam, Vprv, lamprv, pv, pq);

%% linear operator for computing the tangent predictor
J = [   J   dF_dlam; 
      dP_dV dP_dlam ];

% J = [ J, dF_dlam; 
%       z([pv; pq; nb+pq; 2*nb+1])'];

%% compute normalized tangent predictor
z = zprv;
s = zeros(npv+2*npq+1, 1);
s(end,1) = sign(direction);
z([pv; pq; nb+pq; 2*nb+1]) = J\s;   %% tangent vector
z = z/norm(z);                      %% normalize tangent predictor
