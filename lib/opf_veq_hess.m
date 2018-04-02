function d2Veq = opf_veq_hess(x, lambda, mpc, idx, mpopt)
%OPF_VEQ_HESS  Evaluates Hessian of voltage magnitude equality constraint.
%   D2VEQ = OPF_VEQ_HESS(X, LAMBDA, MPC, IDX, MPOPT)
%
%   Hessian evaluation function for voltage magnitudes.
%
%   Inputs:
%     X : optimization vector
%     LAMBDA : column vector of Lagrange multipliers on active and reactive
%              power balance constraints
%     MPC : MATPOWER case struct
%     IDX : index of buses whose voltage magnitudes should be fixed
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     D2VEQ : Hessian of voltage magnitudes.
%
%   Example:
%       d2Veq = opf_veq_hess(x, lambda, mpc, idx, mpopt);
%
%   See also OPF_VEQ_FCN.

%   MATPOWER
%   Copyright (c) 2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% unpack data
[Vr, Vi] = deal(x{:});

%% problem dimensions
nb = length(Vi);            %% number of buses
n = length(idx);            %% number of buses with fixed voltage magnitudes

%% compute voltage magnitude cubed
Vr2 = Vr(idx).^2;
Vi2 = Vi(idx).^2;
VrVi = Vr(idx) .* Vi(idx);
Vm3 = (Vr2 + Vi2).^(3/2);   %% Vm.^3;

%%----- evaluate Hessian of voltage magnitude constraints -----
lamVm_over_Vm3 = lambda ./ Vm3;

Vm_rr = sparse(idx, idx,  Vi2  .* lamVm_over_Vm3, nb, nb);
Vm_ri = sparse(idx, idx, -VrVi .* lamVm_over_Vm3, nb, nb);
Vm_ir = Vm_ri;
Vm_ii = sparse(idx, idx,  Vr2  .* lamVm_over_Vm3, nb, nb);

%% construct Hessian
d2Veq = [Vm_rr Vm_ri; Vm_ir Vm_ii];
