function d2G = opf_current_balance_hess(x, lambda, mpc, Ybus, mpopt)
%OPF_CURRENT_BALANCE_HESS  Evaluates Hessian of current balance constraints.
%   D2G = OPF_CURRENT_BALANCE_HESS(X, LAMBDA, OM, YBUS, MPOPT)
%
%   Hessian evaluation function for AC real and imaginary current balance
%   constraints.
%
%   Inputs:
%     X : optimization vector
%     LAMBDA : column vector of Lagrange multipliers on active and reactive
%              power balance constraints
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
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
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
[Var, Vmi, Pg, Qg] = deal(x{:});

%% problem dimensions
nb = length(Var);           %% number of buses
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

%%----- evaluate Hessian of power balance constraints -----
if mpopt.opf.v_cartesian
    %% reconstruct V
    V = Var + 1j* Vmi;

    %% compute 2nd derivatives
    [Grrr, Grri, Grir, Grii] = d2Imis_dV2_C(Sbus, Ybus, V, lamP);
    [Girr, Giri, Giir, Giii] = d2Imis_dV2_C(Sbus, Ybus, V, lamQ);
    [Gr_SV, Gr_VS] = d2Imis_dVdSg_C(V, lamP, Cg);
    [Gi_SV, Gi_VS] = d2Imis_dVdSg_C(V, lamQ, Cg);
    % RDZ: The above do not need to return both since the 2nd return arg is
    % always the transpose of the first. Just use the transpose directly below.

    %% construct Hessian
    d2G = [
        real([Grrr Grri; Grir Grii]) + imag([Girr Giri; Giir Giii]), real(Gr_VS) + imag(Gi_VS);
        real(Gr_SV) + imag(Gi_SV), sparse(2*ng, 2*ng)
    ];
else
    %% reconstruct V
    V = Vmi .* exp(1j * Var);

    %% compute 2nd derivatives
    [Graa, Grav, Grva, Grvv] = d2Imis_dV2_P(Sbus, Ybus, V, lamP);
    [Giaa, Giav, Giva, Givv] = d2Imis_dV2_P(Sbus, Ybus, V, lamQ);
    [Gr_SV, Gr_VS] = d2Imis_dVdSg_P(V, lamP, Cg);
    [Gi_SV, Gi_VS] = d2Imis_dVdSg_P(V, lamQ, Cg);

    %% construct Hessian
    d2G = [
        real([Graa Grav; Grva Grvv]) + imag([Giaa Giav; Giva Givv]), real(Gr_VS) + imag(Gi_VS);
        real(Gr_SV) + imag(Gi_SV), sparse(2*ng, 2*ng)
    ];
end
