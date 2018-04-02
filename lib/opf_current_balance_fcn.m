function [g, dg] = opf_current_balance_fcn(x, mpc, Ybus, mpopt)
%OPF_CURRENT_BALANCE_FCN  Evaluates AC current balance constraints and their gradients.
%   [G, DG] = OPF_CURRENT_BALANCE_FCN(X, OM, YBUS, MPOPT)
%
%   Computes the real or imaginary current balance equality constraints for
%   AC optimal power flow. Computes constraint vectors and their gradients.
%
%   Inputs:
%     X : optimization vector
%     MPC : MATPOWER case struct
%     YBUS : bus admittance matrix
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     G  : vector of equality constraint values (real/imaginary current balances)
%     DG : (optional) equality constraint gradients
%
%   Examples:
%       g = opf_current_balance_fcn(x, mpc, Ybus, mpopt);
%       [g, dg] = opf_current_balance_fcn(x, mpc, Ybus, mpopt);
%
%   See also OPF_POWER_BALANCE_HESS

%   MATPOWER
%   Copyright (c) 1996-2018, Power Systems Engineering Research Center (PSERC)
%   by Baljinnyam Sereeter, Delft University of Technology
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialize -----
%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% unpack data
[baseMVA, bus, gen] = deal(mpc.baseMVA, mpc.bus, mpc.gen);
if mpopt.opf.v_cartesian
    [Vr, Vi, Pg, Qg] = deal(x{:});
    V = Vr + 1j* Vi;        %% reconstruct V
else
    [Va, Vm, Pg, Qg] = deal(x{:});
    V = Vm .* exp(1j * Va); %% reconstruct V
end

%% problem dimensions
nb = length(V);             %% number of buses
ng = length(Pg);            %% number of dispatchable injections

%% ----- evaluate constraints -----
%% put Pg, Qg back in gen
gen(:, PG) = Pg * baseMVA;  %% active generation in MW
gen(:, QG) = Qg * baseMVA;  %% reactive generation in MVAr

%% rebuild Sbus
Sbus = makeSbus(baseMVA, bus, gen);  %% net injected power in p.u.

%% evaluate complex current balance mismatches
mis = Ybus*V - conj(Sbus./V);

%% assemble real and imaginary current balance constraints
g = [ real(mis);    %% real current mismatch
      imag(mis) ];  %% imaginary current mismatch

%%----- evaluate constraint gradients -----
if nargout > 1
    %% compute partials of injected bus powers
    Cg = sparse(gen(:, GEN_BUS), 1:ng, 1, nb, ng);
    InvConjV = sparse(1:nb, 1:nb, 1./conj(V), nb, nb);
    dImis_dPg = -InvConjV * Cg;     %% dImis w.r.t. Pg
    dImis_dQg = -1j * dImis_dPg;    %% dImis w.r.t. Qg
    [dImis_dV1, dImis_dV2] = dImis_dV(Sbus, Ybus, V, mpopt.opf.v_cartesian);     %% w.r.t. V
    dg = [
        real([dImis_dV1 dImis_dV2 dImis_dPg dImis_dQg]);    %% Ir mismatch w.r.t V1, V2, Pg, Qg
        imag([dImis_dV1 dImis_dV2 dImis_dPg dImis_dQg]) ;   %% Ii mismatch w.r.t V1, V2, Pg, Qg
    ];
end
