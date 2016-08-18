function [V0, lam0] = cpf_predictor(V, lam, z, step, pv, pq)
%CPF_PREDICTOR  Performs the predictor step for the continuation power flow
%   [V0, LAM0] = CPF_PREDICTOR(V, LAM, Z, STEP, PV, PQ)
%
%   Computes a prediction (approximation) to the next solution of the
%   continuation power flow using a normalized tangent predictor.
%
%   Inputs:
%       V : complex bus voltage vector at current solution
%       LAM : scalar lambda value at current solution
%       Z : normalized tangent prediction vector from previous step
%       STEP : continuation step length
%       PV : vector of indices of PV buses
%       PQ : vector of indices of PQ buses
%
%   Outputs:
%       V0 : predicted complex bus voltage vector
%       LAM0 : predicted lambda continuation parameter

%   MATPOWER
%   Copyright (c) 1996-2016 by Power System Engineering Research Center (PSERC)
%   by Shrirang Abhyankar, Argonne National Laboratory
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% sizes
nb = length(V);

Va = angle(V);
Vm = abs(V);
Va0 = Va;
Vm0 = Vm;

%% prediction for next step
Va0([pv; pq]) = Va([pv; pq]) + step * z([pv; pq]);
Vm0([pq])     = Vm([pq])     + step * z([nb+pq]);
lam0 = lam + step * z(2*nb+1);
V0 = Vm0 .* exp(1j * Va0);
