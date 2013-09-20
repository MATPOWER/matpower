function [V0, lam0, z] = cpf_predictor(Vprv, lamprv, Ybus, Sxfr, pv, pq, step, z)
%CPF_PREDICTOR  Performs the predictor step for the continuation power flow
%           by computing normalized tangent predictor and the predicted
%           solution for the next continuation step.
%
% [V0, lam0, z] = cpf_predictor(Vprv, lamprv, Ybus, Sxfr, pv, pq, step, z)
% returns the predicted complex voltage vector V0, predicted loading
% parameter lam0 for the next continuation step. cpf_predictor() also returns
% the normalized tangent predictor z.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2013 by Power System Engineering Research Center (PSERC)
%
%   Contributed by Shrirang Abhyankar

j = sqrt(-1);
nb = length(Vprv);
npv = length(pv);
npq = length(pq);

%% Compute Jacobian for the power flow equations
[dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, Vprv);
    
j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
j12 = real(dSbus_dVm([pv; pq], pq));
j21 = imag(dSbus_dVa(pq, [pv; pq]));
j22 = imag(dSbus_dVm(pq, pq));

J = [   j11 j12;
        j21 j22;    ];

%% Linear operator for computing the tangent predictor
J = [ J, -[real(Sxfr([pv; pq])); imag(Sxfr(pq))]; 
      z([pv; pq; nb+pq; 2*nb+1])'];

Vaprv=angle(Vprv);
Vmprv=abs(Vprv);

%% Compute normalized tangent predictor
s=zeros(npv+2*npq+1,1);
s(end,1)= 1;                     %% Increase in the direction of lambda
z([pv; pq; nb+pq; 2*nb+1])= J\s; %% tangent vector
z = z/norm(z);    % normalize tangent predictor

Va0 = Vaprv;
Vm0 = Vmprv;
lam0 = lamprv;

%% prediction for next step
Va0([pv; pq])=Vaprv([pv; pq])+ step* z([pv; pq]); 
Vm0([pq])=Vmprv([pq]) +step* z([nb+pq]);
lam0 = lamprv +step*z(2*nb+1);
V0 = Vm0.*exp(j*Va0);