function [g, dg] = opf_power_balance_fcn(x, mpc, Ybus, mpopt)
%OPF_POWER_BALANCE_FCN  Evaluates AC power balance constraints and their gradients.
%   [G, DG] = OPF_POWER_BALANCE_FCN(X, OM, YBUS, MPOPT)
%
%   Computes the active or reactive power balance equality constraints for
%   AC optimal power flow. Computes constraint vectors and their gradients.
%
%   Inputs:
%     X : optimization vector
%     MPC : MATPOWER case struct
%     YBUS : bus admittance matrix
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     G  : vector of equality constraint values (active/reactive power balances)
%     DG : (optional) equality constraint gradients
%
%   Examples:
%       g = opf_power_balance_fcn(x, mpc, Ybus, mpopt);
%       [g, dg] = opf_power_balance_fcn(x, mpc, Ybus, mpopt);
%
%   See also OPF_POWER_BALANCE_HESS

%   MATPOWER
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
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

%% ----- evaluate constraints -----
%% put Pg, Qg back in gen
gen(:, PG) = Pg * baseMVA;  %% active generation in MW
gen(:, QG) = Qg * baseMVA;  %% reactive generation in MVAr

%% rebuild Sbus
if mpopt.opf.v_cartesian
    Sbus = makeSbus(baseMVA, bus, gen);             %% net injected power in p.u.
else
    Sbus = makeSbus(baseMVA, bus, gen, mpopt, Vm);  %% net injected power in p.u.
end

%% evaluate complex power balance mismatches
mis = V .* conj(Ybus * V) - Sbus;

%% assemble active and reactive power balance constraints
g = [ real(mis);    %% active power mismatch
      imag(mis) ];  %% reactive power mismatch

%%----- evaluate constraint gradients -----
if nargout > 1
    %% compute partials of injected bus powers
    [dSbus_dV1, dSbus_dV2] = dSbus_dV(Ybus, V, mpopt.opf.v_cartesian);  %% w.r.t. V
    neg_Cg = sparse(gen(:, GEN_BUS), 1:ng, -1, nb, ng);     %% Pbus w.r.t. Pg
    if ~mpopt.opf.v_cartesian
        %% adjust for voltage dependent loads
        [dummy, neg_dSd_dVm] = makeSbus(baseMVA, bus, gen, mpopt, Vm);
        dSbus_dV2 = dSbus_dV2 - neg_dSd_dVm;
    end
    dg = [
        real([dSbus_dV1 dSbus_dV2]) neg_Cg sparse(nb, ng);  %% P mismatch w.r.t V1, V2, Pg, Qg
        imag([dSbus_dV1 dSbus_dV2]) sparse(nb, ng) neg_Cg;  %% Q mismatch w.r.t V1, V2, Pg, Qg
    ];
end
