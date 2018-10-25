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

%%----- evaluate Hessian of power balance constraints -----
%% compute 2nd derivatives
[Gp11, Gp12, Gp21, Gp22] = d2Sbus_dV2(Ybus, V, lamP, mpopt.opf.v_cartesian);
[Gq11, Gq12, Gq21, Gq22] = d2Sbus_dV2(Ybus, V, lamQ, mpopt.opf.v_cartesian);

if ~mpopt.opf.v_cartesian
    %% adjust for voltage dependent loads (constant impedance part of ZIP loads)
    diaglam = sparse(1:nb, 1:nb, lamP, nb, nb);
    Sd = makeSdzip(mpc.baseMVA, mpc.bus, mpopt);
    diagSdz = sparse(1:nb, 1:nb, Sd.z, nb, nb);
    Gp22 = Gp22 + 2 * diaglam * diagSdz;
end

%% construct Hessian
d2G = [
    real([Gp11 Gp12; Gp21 Gp22]) + imag([Gq11 Gq12; Gq21 Gq22]) sparse(2*nb, 2*ng);
    sparse(2*ng, 2*nb + 2*ng)
];
