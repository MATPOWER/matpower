function d2G = opf_current_balance_hess(x, lambda, mpc, Ybus, mpopt)
%OPF_CURRENT_BALANCE_HESS  Evaluates Hessian of current balance constraints.
%   D2G = OPF_CURRENT_BALANCE_HESS(X, LAMBDA, OM, YBUS, MPOPT)
%
%   Hessian evaluation function for AC real and imaginary current balance
%   constraints.
%
%   Inputs:
%     X : optimization vector
%     LAMBDA : column vector of Lagrange multipliers on real and imaginary
%              current balance constraints
%     MPC : MATPOWER case struct
%     YBUS : bus admittance matrix
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     D2G : Hessian of current balance constraints.
%
%   Example:
%       d2G = opf_current_balance_hess(x, lambda, mpc, Ybus, mpopt);
%
%   See also OPF_CURRENT_BALANCE_FCN.

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
    V = Vr + 1j * Vi;           %% reconstruct V
else
    [Va, Vm, Pg, Qg] = deal(x{:});
    V = Vm .* exp(1j * Va);     %% reconstruct V
end

%% problem dimensions
nb = length(V);             %% number of buses
ng = length(Pg);            %% number of dispatchable injections

nlam = length(lambda) / 2;
lamP = lambda(1:nlam);
lamQ = lambda((1:nlam)+nlam);
Cg = sparse(gen(:, GEN_BUS), 1:ng, 1, nb, ng);

%% put Pg, Qg back in gen
gen(:, PG) = Pg * baseMVA;  %% active generation in MW
gen(:, QG) = Qg * baseMVA;  %% reactive generation in MVAr

%% rebuild Sbus
Sbus = makeSbus(baseMVA, bus, gen);

%%----- evaluate Hessian of current balance constraints -----
%% compute 2nd derivatives
[Gr11, Gr12, Gr21, Gr22] = d2Imis_dV2(Sbus, Ybus, V, lamP, mpopt.opf.v_cartesian);
[Gi11, Gi12, Gi21, Gi22] = d2Imis_dV2(Sbus, Ybus, V, lamQ, mpopt.opf.v_cartesian);
Gr_sv = d2Imis_dVdSg(Cg, V, lamP, mpopt.opf.v_cartesian);
Gi_sv = d2Imis_dVdSg(Cg, V, lamQ, mpopt.opf.v_cartesian);

%% construct Hessian
d2G = [
    real([Gr11 Gr12; Gr21 Gr22]) + imag([Gi11 Gi12; Gi21 Gi22]), real(Gr_sv.') + imag(Gi_sv.');
    real(Gr_sv) + imag(Gi_sv), sparse(2*ng, 2*ng)
];
