function d2G = opf_power_balance_hess(x, lambda, mpc, Ybus, mpopt)
%OPF_POWER_BALANCE_HESS  Evaluates Hessian of power balance constraints.
%   D2G = OPF_POWER_BALANCE_HESS(X, LAMBDA, OM, YBUS, MPOPT)
%
%   Hessian evaluation function for AC active and reactive power balance
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
%     D2G : Hessian of power balance constraints.
%
%   Example:
%       d2G = opf_power_balance_hess(x, lambda, mpc, Ybus, mpopt);
%
%   See also OPF_POWER_BALANCE_FCN.

%   MATPOWER
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialize -----
%% unpack data
[Var, Vmi, Pg, Qg] = deal(x{:});

%% problem dimensions
nb = length(Var);           %% number of buses
ng = length(Pg);            %% number of dispatchable injections

nlam = length(lambda) / 2;
lamP = lambda(1:nlam);
lamQ = lambda((1:nlam)+nlam);

%%----- evaluate Hessian of power balance constraints -----
if mpopt.opf.v_cartesian
    %% reconstruct V
    V = Var + 1j* Vmi;

    %% compute 2nd derivatives
    [Gprr, Gpri, Gpir, Gpii] = d2Sbus_dV2_C(Ybus, V, lamP);
    [Gqrr, Gqri, Gqir, Gqii] = d2Sbus_dV2_C(Ybus, V, lamQ);

    %% construct Hessian
    d2G = [
        real([Gprr Gpri; Gpir Gpii]) + imag([Gqrr Gqri; Gqir Gqii]) sparse(2*nb, 2*ng);
        sparse(2*ng, 2*nb + 2*ng)
    ];
else
    %% reconstruct V
    V = Vmi .* exp(1j * Var);

    %% compute 2nd derivatives
    [Gpaa, Gpav, Gpva, Gpvv] = d2Sbus_dV2_P(Ybus, V, lamP);
    [Gqaa, Gqav, Gqva, Gqvv] = d2Sbus_dV2_P(Ybus, V, lamQ);

    %% adjust for voltage dependent loads (constant impedance part of ZIP loads)
    diaglam = sparse(1:nb, 1:nb, lamP, nb, nb);
    Sd = makeSdzip(mpc.baseMVA, mpc.bus, mpopt);
    diagSdz = sparse(1:nb, 1:nb, Sd.z, nb, nb);
    Gpvv = Gpvv + 2 * diaglam * diagSdz;

    %% construct Hessian
    d2G = [
        real([Gpaa Gpav; Gpva Gpvv]) + imag([Gqaa Gqav; Gqva Gqvv]) sparse(2*nb, 2*ng);
        sparse(2*ng, 2*nb + 2*ng)
    ];
end
