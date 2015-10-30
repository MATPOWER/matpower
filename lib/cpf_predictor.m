function [V0, lam0, z] = cpf_predictor(V, lam, Ybus, Sxfr, Sbust, pv, pq, ...
                            step, z, Vprv, lamprv, parameterization)
%CPF_PREDICTOR  Performs the predictor step for the continuation power flow
%   [V0, LAM0, Z] = CPF_PREDICTOR(V, LAM, YBUS, SXFR, PV, PQ, STEP, Z, ...
%                                 VPRV, LAMPRV, PARAMETERIZATION)
%
%   Computes a prediction (approximation) to the next solution of the
%   continuation power flow using a normalized tangent predictor.
%
%   Inputs:
%       V : complex bus voltage vector at current solution
%       LAM : scalar lambda value at current solution
%       YBUS : complex bus admittance matrix
%       SXFR : handle of function returning complex vector of scheduled
%              transfers in p.u. (difference between bus injections in base
%              and target cases)
%       SBUST : handle of function returning bus injections for target case
%               and derivatives w.r.t. voltage magnitudes
%       PV : vector of indices of PV buses
%       PQ : vector of indices of PQ buses
%       STEP : continuation step length
%       Z : normalized tangent prediction vector from previous step
%       VPRV : complex bus voltage vector at previous solution
%       LAMPRV : scalar lambda value at previous solution
%       PARAMETERIZATION : Value of cpf.parameterization option.
%
%   Outputs:
%       V0 : predicted complex bus voltage vector
%       LAM0 : predicted lambda continuation parameter
%       Z : the normalized tangent prediction vector

%   MATPOWER
%   Copyright (c) 1996-2015 by Power System Engineering Research Center (PSERC)
%   by Shrirang Abhyankar, Argonne National Laboratory
%   and Ray Zimmerman, PSERC Cornell
%
%   $Id$
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
[dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V);
[dummy, neg_dSd_dVm] = Sbust(Vm);
dSbus_dVm = dSbus_dVm - lam * neg_dSd_dVm;

j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
j12 = real(dSbus_dVm([pv; pq], pq));
j21 = imag(dSbus_dVa(pq, [pv; pq]));
j22 = imag(dSbus_dVm(pq, pq));

J = [   j11 j12;
        j21 j22;    ];

Sxf = Sxfr(Vm);     %% update target based on current voltages
dF_dlam = -[real(Sxf([pv; pq])); imag(Sxf(pq))];
[dP_dV, dP_dlam] = cpf_p_jac(parameterization, z, V, lam, Vprv, lamprv, pv, pq);

%% linear operator for computing the tangent predictor
J = [   J   dF_dlam; 
      dP_dV dP_dlam ];

% J = [ J, dF_dlam; 
%       z([pv; pq; nb+pq; 2*nb+1])'];

Vaprv = angle(V);
Vmprv = abs(V);

%% compute normalized tangent predictor
s = zeros(npv+2*npq+1, 1);
s(end,1) = 1;                       %% increase in the direction of lambda
z([pv; pq; nb+pq; 2*nb+1]) = J\s;   %% tangent vector
z = z/norm(z);                      %% normalize tangent predictor

Va0 = Vaprv;
Vm0 = Vmprv;
lam0 = lam;

%% prediction for next step
Va0([pv; pq]) = Vaprv([pv; pq]) + step * z([pv; pq]);
Vm0([pq]) = Vmprv([pq]) + step * z([nb+pq]);
lam0 = lam + step * z(2*nb+1);
V0 = Vm0 .* exp(1j * Va0);
