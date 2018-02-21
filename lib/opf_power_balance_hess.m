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
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% ----- evaluate Hessian of power balance constraints -----
    nlam = length(lambda) / 2;
    lamP = lambda(1:nlam);
    lamQ = lambda((1:nlam)+nlam);
    
%%----- initialize -----    
if mpopt.opf.v_cartesian
    [Vi, Vr, Pg, Qg] = deal(x{:});  %% unpack data
    %% reconstruct V
    V = Vr + 1j* Vi;  
    %% problem dimensions
    nb = length(Vi);            %% number of buses
    ng = length(Pg);            %% number of dispatchable injections
    
    [Gpii, Gpir, Gpri, Gprr] = d2Sbus_dV2_C(Ybus, V, lamP);
    [Gqii, Gqir, Gqri, Gqrr] = d2Sbus_dV2_C(Ybus, V, lamQ);    
    %% construct Hessian
    d2G = [
        real([Gpii Gpir; Gpri Gprr]) + imag([Gqii Gqir; Gqri Gqrr]) sparse(2*nb, 2*ng);
        sparse(2*ng, 2*nb + 2*ng)
    ];      
else
    [Va, Vm, Pg, Qg] = deal(x{:});  %% unpack data
    %% reconstruct V
    V = Vm .* exp(1j * Va);
    
    %% problem dimensions
    nb = length(Va);            %% number of buses
    ng = length(Pg);            %% number of dispatchable injections

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


